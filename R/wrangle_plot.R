#wrangle for plot

# -----------------------------------------------------------------------------
# format_to_id_map()
# Why:
# - downscalR and your tables sometimes use "ns" as a character id (id_c).
# - For raster painting you use an integer id (ns_int) stored in raster cells.
# This helper converts res$out.res$ns (character) -> integer ns_int consistently.
# -----------------------------------------------------------------------------
format_to_id_map <- function(res, ns_map) {
  
  temp <- res$out.res
  
  formatted <- temp %>%
    dplyr::mutate(ns = as.character(ns)) %>%
    dplyr::left_join(ns_map %>% dplyr::select(ns, ns_int), by = "ns") %>%
    dplyr::mutate(ns = ns_int) %>%
    dplyr::select(-ns_int)
  
  formatted
}

# -----------------------------------------------------------------------------
# loss_and_gains()
# Why:
# - downscalR results are stored as transition flows: (lu.from -> lu.to) per cell and year.
# - For “change maps” you usually want *net gains/losses per final land-use class*,
#   without caring which specific origin class was replaced.
#
# What it returns:
# - A long table with one row per (ns, times, lu.to) giving:
#   value > 0  => net gain of that land-use in that cell and year
#   value < 0  => net loss of that land-use in that cell and year
#
# How:
# 1) Gains: sum all transitions ending in each lu.to (excluding diagonal lu.from==lu.to)
# 2) Losses: sum all transitions leaving each lu.from, then multiply by -1
# 3) Combine and collapse to a single signed value per (lu, ns, times)
# -----------------------------------------------------------------------------
loss_and_gains <- function(res) {
  
  temp <- res$out.res
  
  to.plot <- temp %>%
    dplyr::filter(lu.from != lu.to) %>%
    # Gains by destination class
    dplyr::group_by(lu.to, ns, times) %>%
    dplyr::summarise(gain = sum(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Type = "gain", lu = lu.to) %>%
    
    # Bind losses by origin class (negative sign for map interpretation)
    dplyr::bind_rows(
      temp %>%
        dplyr::filter(lu.from != lu.to) %>%
        dplyr::group_by(lu.from, ns, times) %>%
        dplyr::summarise(loss = -sum(value, na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(Type = "loss", lu = lu.from)
    ) %>%
    
    # Collapse gains and losses into one signed number per cell and land-use
    dplyr::group_by(lu, times, ns) %>%
    dplyr::summarise(value = sum(gain, loss, na.rm = TRUE), .groups = "drop") %>%
    dplyr::rename(lu.to = lu)
  
  return(to.plot)
}


# -----------------------------------------------------------------------------
# cum_difference()
# Why:
# - downscalR outputs can be interpreted as land-use “levels” per year
#   (area in each land-use class per cell at that time).
# - A common diagnostic is “how much did each land-use change between 2020 and 2050?”
#   per cell and final land-use class.
#
# What it returns:
# - A long table with one row per (ns, lu.to) giving:
#   value = area_2050 - area_2020
#   times is set to 2050 just to keep a single “time” label for plotting functions.
#
# How:
# 1) Keep only years 2020 and 2050
# 2) Aggregate within each (ns, lu.to, times) to get total area by land-use
# 3) Reshape wide to get columns `2020` and `2050`
# 4) Compute the difference and store it as value
# -----------------------------------------------------------------------------
cum_difference <- function(res) {
  
  temp <- res$out.res
  
  to.plot <- temp %>%
    dplyr::filter(times %in% c(2020, 2050)) %>%
    dplyr::group_by(lu.to, ns, times) %>%
    dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = times, values_from = value) %>%
    dplyr::mutate(value = `2050` - `2020`,
                  times = 2050)
  
  return(to.plot)
}



# -----------------------------------------------------------------------------
# pairwise_net_difference()
# Why:
# - “Net gain/loss by land-use” hides *which land-use replaced which*.
# - Sometimes you want a transition-focused diagnostic: for each pair (A,B),
#   what is the *net* flow between them over the whole horizon?
#
# What it returns:
# - A table with one row per (ns, lu.from, lu.to) representing the *net* transition
#   between that unordered pair, expressed as a percentage of cell area.
#
# Key definition (net between A and B):
# - net(A,B) = flow(A -> B) - flow(B -> A)
# - If net > 0: overall movement was from A to B
# - If net < 0: overall movement was from B to A
#
# How:
# 1) Drop diagonal transitions (lu.from==lu.to)
# 2) For each row, define an ordered pair (a,b) where a=min(from,to), b=max(from,to)
# 3) Assign sign = +1 if the flow is a->b, -1 if the flow is b->a
# 4) Sum sign*value across all years to get the net flow for that (ns,a,b)
# 5) Convert that net flow into % of pixel area using the grid area
#
# -----------------------------------------------------------------------------
pairwise_net_difference <- function(res, ns_map, grid){
  
  temp <- res$out.res
  
  to.plot <- temp %>%
    #do not keep 
    filter(lu.from != lu.to) %>%
    #Positive = net flow from the alphabetically-first class (a) to the second (b); 
    #negative = net the other way.
    #will be used to compute the net flow from a to b
    mutate(
      a = pmin(lu.from, lu.to),
      b = pmax(lu.from, lu.to),
      sign = if_else(lu.from == a & lu.to == b, 1, -1)
    ) %>% 
    #adding the missing b to a flow 
    rbind(
      temp %>%
        filter(lu.from != lu.to) %>%
        mutate(
          b = pmin(lu.from, lu.to),
          a = pmax(lu.from, lu.to),
          sign = if_else(lu.from == a & lu.to == b, 1, -1)
        )
    ) %>%
    #compute the net flow for each pair
    group_by(ns, a, b) %>%
    summarise(net = sum(sign * value, na.rm = TRUE), .groups = "drop") %>%
    arrange( ns, a, b) %>% 
    rename(lu.from = a, 
           lu.to   = b) %>% 
    left_join(ns_map %>% select(-ns), by = c("ns" = "ns_int")) %>% 
    left_join(grid %>% select(id_c, area)) %>% 
    mutate(value = 100 * net/(area/1000)) %>% 
    select(ns, lu.from, lu.to, value) %>% 
    filter(lu.from != "urban") %>% 
    filter(lu.to != "forest") %>% 
    filter(lu.from != "newforest") 
  
  return(to.plot)
}

pairwise_net_difference <- function(res, ns_map) {
  
  temp <- res$out.res
  
  to.plot <- temp %>%
    dplyr::filter(lu.from != lu.to) %>%
    
    # Create a consistent ordering of each pair so A<->B are grouped together
    dplyr::mutate(
      a    = pmin(lu.from, lu.to),
      b    = pmax(lu.from, lu.to),
      sign = dplyr::if_else(lu.from == a & lu.to == b, 1, -1)
    ) %>%
    
    # Duplicate rows with swapped labels to ensure both directions
    # are represented consistently before aggregation.
    dplyr::bind_rows(
      temp %>%
        dplyr::filter(lu.from != lu.to) %>%
        dplyr::mutate(
          b    = pmin(lu.from, lu.to),
          a    = pmax(lu.from, lu.to),
          sign = dplyr::if_else(lu.from == a & lu.to == b, 1, -1)
        )
    ) %>%
    
    # Net flow for each unordered pair (a,b) in each cell
    dplyr::group_by(ns, a, b) %>%
    dplyr::summarise(net = sum(sign * value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(ns, a, b) %>%
    dplyr::rename(lu.from = a, lu.to = b) %>%
    
    # Attach original id_c (optional) and cell area for % conversion
    dplyr::left_join(ns_map %>% dplyr::select(-ns), by = c("ns" = "ns_int")) %>%
    dplyr::left_join(grid %>% dplyr::select(id_c, area), by = "id_c") %>%
    
    # Convert to % of pixel:
    # area is in km2 in your grid; (area/1000) is thousand ha
    # multiplying by 100 expresses "percentage of pixel"
    dplyr::mutate(value = 100 * net / (area / 1000)) %>%
    dplyr::select(ns, lu.from, lu.to, value) %>%
    
    # Project-specific filters (keep/remove depending on your reporting needs)
    dplyr::filter(lu.from != "urban") %>%
    dplyr::filter(lu.to   != "forest") %>%
    dplyr::filter(lu.from != "newforest")
  
  return(to.plot)
}
