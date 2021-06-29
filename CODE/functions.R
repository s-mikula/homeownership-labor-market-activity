#### House Coordinates ####
# Takes IDOB (house ID) defined as character or numeric and returns house coordinates in S-JTSK.
# Data are retrieved from RSO database at CSO website (czso.cz)
#
# Arguments:
# IDOB ... house ID
#
# Values
# a tibble with columns IDOB, SX, and SY

get_XY <- function(IDOB){
  require(dplyr)
  require(rvest)
  require(xml2)
  require(stringr)
  
  # You need a house ID used internally in the CSO website first
  
  # A) Get an url of the website with general information on the house defined by IDOB
  query_rso <- str_c(
    "http://apl.czso.cz/irso4/budlist.jsp?b=12&idob=",
    IDOB,
    "&pcbudov=&ruianso_id=&ruiantea_id=&cuzkbud_id=&idadr=&ruianam_id=&pcd=&hledej=Vyhledat"
  )
  
  # B) Get a website and isolate link to the map
  read_html(query_rso) %>% 
    # Get table
    html_node(xpath = '//*[@id="core"]/div[3]/table') %>% 
    html_nodes("a") %>% 
    as.character() %>% 
    str_subset("onclick") %>% 
    str_split(" ") %>% 
    str_subset("onclick") -> ids_string
  
  # C) If there are multiple hits (with identical IDOB) use the first one
  if(length(ids_string)>1){
    ids_string[1] -> ids_string
    warning(str_c("IDOB ",IDOB," má více záznamů. Rohový dům?"))
  }
  
  # D) Extract internal IDs
  budid <- ids_string %>% str_extract("budid=(\\d*)") %>% str_replace_all("\\D*","")
  obrprvid <- ids_string %>% str_extract("obrprvid=(\\d*)") %>% str_replace_all("\\D*","")
  
  # Retrieve website with map and marked location of the house
  
  # A) Create an url
  query_map <- str_c(
    "http://apl.czso.cz/irso4/mapa.jsp?pid=1&budid=",
    budid,
    "&obrprvid=",
    obrprvid
  )
  
  # B) Get the website and isolate coordinates
  read_html(query_map) %>% 
    html_node(xpath = '/html/head/script[7]') %>% 
    as.character() %>% 
    str_extract("init\\(-(\\d*),-(\\d*)") %>% 
    str_replace("init\\(","") %>% 
    str_split(",") %>% unlist %>% 
    str_replace_all(" ","") %>% 
    as.integer() -> ids_integer
  
  # Create the output tibble
  
  data_frame(
    IDOB = as.character(IDOB),
    SX = ids_integer[1],
    SY = ids_integer[2]
  )
}

#### Convert coordinates from S-JTSK into WGS84 ####
# Converts S-JTSK into WGS84 using GDAL library

# Arguments
# data ... a data.frame
# SX ... a column with X-coordinates
# SY ... a column with Y-coordinates

# Values
# the input data.frame with additional columns with longitude (lon)
# and latitude (lat)

jtsk2wgs <- function(data, SX = "SX", SY = "SY"){
  require(rgdal)
  require(dplyr)
  
  # Just a quick check
  if(SX == "lon"){
    data <- data %>% rename(lon_orig = lon)
    SX <- "lon_orig"
  }
  
  if(SY == "lat"){
    data <- data %>% rename(lat_orig = lat)
    SY <- "lat_orig"
  }
  
  # Create backup data.frame
  x <- data
  
  coordinates(data) <- c(SX,SY)
  proj4string(data) <- CRS("+init=epsg:5514")
  
  CRS_new <- CRS("+init=epsg:4326")
  data <- spTransform(data, CRS_new)
  
  data@coords %>% as.data.frame() %>% 
    rename(
      lon = SX,
      lat = SY
    ) -> x_data
  
  if(nrow(x_data)!=nrow(x)){
    warning("Size of input and output table differ. Returning transformed S4 spatial object.")
    return(data)
  }else{
    bind_cols(x,x_data)
  }
}

#### Matching ####
#### Wrapper for genetic matching from Matching package -- function is derived from match_towns() 
# from sudetenland project

match_houses <- function(match_dataset, MATCH = NULL, defined_treatment = NULL,
                         caliper = 0.1, replace = FALSE,
                         pop_size = 100, CPUs = "all"){
  # Required libraries
  require(Matching)
  
  if(is.null(MATCH)) MATCH <- c("SX","SY","DUMPOCBYT")
  
  if(is.null(defined_treatment)){
    match_dataset %>% 
      dplyr::mutate(
        treatment = status=="privatization"
      ) -> match_dataset
  }
  
  # Keep only relevant variables
  match_dataset %>% 
    dplyr::select(LIDAMSIDOB,treatment,one_of(MATCH)) ->  dataset
  
  # Keep only complete rows
  dataset <- dataset[complete.cases(dataset),]
  
  
  # Add actual number of cores available
  if(CPUs == "all"){
    require(parallel)
    CPUs <- parallel::detectCores()
  }
  
  # Start cluster
  if(CPUs > 1){
    library(snow)
    cl <- makeCluster(rep("localhost",CPUs), type = "SOCK")
  }
  
  # Get weights
  if(CPUs==1){
    GenMatch(
      Tr = dataset$treatment,
      X = dataset[,MATCH],
      pop.size = pop_size,
      caliper = caliper,
      replace = replace
    ) -> W
  }else{
    try(
      GenMatch(
        Tr = dataset$treatment,
        X = dataset[,MATCH],
        pop.size = pop_size,
        caliper = caliper,
        replace = replace,
        cluster = cl 
      ) -> W
    )
  }
  
  if(CPUs > 1) stopCluster(cl)
  if(!exists("W")) stop("Problem with GenMatch(). Clusters probably?") 
  
  # Get matches
  Match(
    Tr = dataset$treatment,
    X = dataset[,MATCH],
    Weight.matrix = W,
    caliper = caliper,
    replace = replace
  ) -> match
  
  
  detach("package:Matching", unload = TRUE)
  detach("package:MASS", unload = TRUE)
  
  dplyr::bind_cols(
    dataset[match$index.treated,"LIDAMSIDOB"] %>% rename(treatment = LIDAMSIDOB),
    dataset[match$index.control,"LIDAMSIDOB"] %>% rename(control = LIDAMSIDOB)
  ) -> matches
  
  list(
    caliper = caliper,
    no_of_matches = nrow(matches),
    variables_matched_on = MATCH,
    data = match_dataset,
    matches = matches,
    output_GenMatch = W,
    output_Match = match
  )
}

#### Sample ####
get_sample_maxdist <- function(max_dist = NULL, ID_control = control, ID_treatment = treatment, distances_table = distances){
  require(dplyr)
  
  if(is.null(max_dist)) stop("maxdist must be set!")
  
  if(max_dist > 0){
  distances_table %>% 
    filter(ID1 %in% ID_treatment, ID2 %in% ID_control) %>% 
    rename(treatment = ID1, control = ID2) %>% 
    filter(distance <= max_dist)
  }else{
  distances_table %>% 
    filter(ID1 %in% ID_treatment, ID2 %in% ID_control) %>% 
    rename(treatment = ID1, control = ID2) %>% 
    group_by(treatment) %>% 
    arrange(distance) %>% 
    slice(1) %>% 
    ungroup 
  }
}

add_to_sample <- function(x, column_name, data = all_samples){
  x %>% 
    select(control,treatment) %>% 
    gather(group,LIDAMSIDOB) %>% 
    distinct() -> x
  names(x)[1] <- column_name
  
  left_join(data,x,by="LIDAMSIDOB")
}

#### Map of Brno ####
# brno_map <- function(){
#   require(dplyr)
#   require(ggplot2)
#   require(rgdal)
#   require(broom)
#   
#   readOGR("DATA/brno_poly.shx") -> brno_poly
#   
#   brno_poly@data %>% as_data_frame() %>% 
#     mutate(
#       id = row_number() - 1
#     ) -> brno_poly_tab
#   
#   brno_poly_tab %>% 
#     filter(name == "Brno" | landuse == "residential") %>% 
#     select(id,osm_id,name) -> brno_poly_tab
#   
#   brno_poly %>% tidy() %>% as_data_frame() %>% rename(lon = long) -> brno_poly
#   
#   brno_poly %>% filter(id %in% brno_poly_tab$id) -> brno_poly
#   
#   readOGR("DATA/brno_routes.shx") %>% tidy() %>% as_data_frame() %>% rename(lon = long) -> brno_lines
#   
#   
#   ggplot(
#     brno_poly,
#     aes(x=lon,y=lat)
#   ) +
#     geom_polygon(
#       data = brno_poly %>% filter(id!=164),
#       aes(group=group),
#       fill = "grey90",
#       alpha = 1
#     ) +
#     geom_path(
#       data = brno_lines,
#       aes(group=group),
#       #color="red",
#       size = 0.4
#     ) +
#     geom_polygon(
#       data = brno_poly %>% filter(id==164),
#       aes(group=group),
#       color = "black",
#       fill = NA,
#       size = 1.5
#     ) +
#     theme_void() +
#     coord_map(xlim = c(16.42,16.73), ylim = c(49.1,49.3)) -> p
#   
#   return(p)
# }
# 
# ### run it and save resulting object
# if(!file.exists("DATA/brno_map_figure.Rdata")){
#   brno_map() -> brno_map_figure
#   save(brno_map_figure, file = "DATA/brno_map_figure.Rdata")
# }

#### Functions from differents projct file ####

get_mm_formula <- function(x,gvar){
  require(stringr)
  str_c(gvar, "~", x, " + 0") %>% as.formula()
}

get_test_formula <- function(x, gvar = "value"){
  str_c(x, " ~ ", gvar) %>% as.formula()
}

get_model_table <- function(mm_formula, data){
  yvar <- model.frame(mm_formula, data = data)[,1] %>% as_tibble()
  xvar <- model.matrix(mm_formula, data = data) %>% as_tibble()
  
  bind_cols(yvar,xvar)
}

get_stats <- function(l){
  #if(is.null(w)){
  l %>% 
    gather(yvariable, yvalues, -1) %>% 
    group_by(value, yvariable) %>% 
    summarise(
      mean = mean(yvalues),
      sd = sd(yvalues),
      n = n()
    ) %>% 
    ungroup()
  # }else{
  #   require(Weighted.Desc.Stat)
  #   
  #   l %>% 
  #     gather(yvariable, yvalues, -1, -one_of(w)) %>% 
  #     group_by(value, yvariable) %>% 
  #     summarise(
  #       mean = weighted.mean(yvalues, w),
  #       sd = w.sd(yvalues, w),
  #       n = n()
  #     ) %>% 
  #     ungroup()
  # }
}

get_stats.0g <- function(l){
  l %>% 
    gather(yvariable, yvalues) %>% 
    group_by(yvariable) %>% 
    summarise(
      mean = mean(yvalues),
      sd = sd(yvalues),
      n = n()
    ) %>% 
    ungroup()
}

get_tests <- function(l){
  require(broom)
  # Get t-test p-value
  lapply(names(l)[-1], get_test_formula) %>% 
    lapply(lm, data = l) %>% 
    lapply(coeftest, vcov. = vcovHC) %>% 
    lapply(tidy) %>% 
    lapply(filter, str_detect(term, "value")) %>% 
    bind_rows() %>% 
    select(p.value) %>% 
    bind_cols(
      data_frame(yvariable = names(l)[-1]),
      .
    )
}

get_tests3 <- function(l){
  require(broom)
  # Get F-test p-value
  lapply(names(l)[-1], get_test_formula) %>% 
    lapply(lm, data = l) %>% 
    lapply(glance) %>% 
    #lapply(filter, str_detect(term, "value")) %>% 
    bind_rows() %>% 
    select(p.value) %>% 
    bind_cols(
      data_frame(yvariable = names(l)[-1]),
      .
    )
}

add_stars <- function(pval){
  y <- character(length(pval))
  y[pval < 0.1]  <- "*  "
  y[pval < 0.05] <- "** "
  y[pval < 0.01] <- "***"
  
  return(y)
}

get_comp_table <- function(xtab, stars = TRUE){
  require(dplyr)
  require(stringr)
  
  # Is the input a data.frame?
  if(!is.data.frame(xtab)){
    stop("Input has to be a data.frame.")
  }
  
  # If there is grouping with only one group -> ungroup()
  if(is.grouped_df(xtab) & length(unique(group_indices(xtab))) == 1){
    message("Table contains only one group. Ungrouping.")
    xtab <- xtab %>% ungroup()
  }
  
  
  if(is.grouped_df(xtab)){
    # Code for grouped tables
    
    # List of grouping variables
    gvars <- group_vars(xtab)
    
    # List of other variables (to be compared)
    list_of_vars <- names(xtab)
    list_of_vars <- list_of_vars[!(list_of_vars %in% gvars)]
    
    # Create tables separataly for each variable
    list_of_vars %>% 
      lapply(get_mm_formula, gvar = gvars) %>% 
      lapply(get_model_table, data = xtab) -> ll
    
    # Create table with summary stats
    ll %>% 
      lapply(get_stats) %>% 
      bind_rows() -> ll_stats
    
    
    if(length(unique(group_indices(xtab))) == 2){
      # t-test on means (2 level grouping)
      ll %>% 
        lapply(get_tests) %>% 
        bind_rows() %>% 
        left_join(ll_stats,., by = "yvariable") %>% 
        arrange(yvariable,value) -> out
    }else{
      # f-test (>2 level grouping)
      ll %>% 
        lapply(get_tests3) %>% 
        bind_rows() %>% 
        left_join(ll_stats,., by = "yvariable") %>% 
        arrange(yvariable,value) -> out
      message("Returning F-test p-values.")
    }
    
    if(stars){
      out %>% 
        mutate(
          stars = add_stars(p.value)
        ) -> out
    }
    
  }else{
    # Code for ungrouped variables
    
    names(xtab) %>% 
      lapply(get_mm_formula, gvar = "") %>% 
      lapply(get_model_table, data = xtab) %>% 
      lapply(select, -value) -> ll
    
    ll %>% 
      lapply(get_stats.0g) %>% 
      bind_rows() -> out
    
    message("No groups defined. Use group_by().")
  }
  
  return(out)
}

find_matches <- function(house_ID, h1 = r91_data, h2 = r01_data){
  h1 <- h1 %>% filter(adresado_01_klic == house_ID)
  h2 <- h2 %>% filter(adresado_01_klic == house_ID)
  
  crossing(
    ID_91 = h1$person_ID_91,
    ID_01 = h2$person_ID
  ) %>% 
    rowwise() %>% 
    mutate(
      sex = h1$POHLAVI[h1$person_ID_91 == ID_91] == h2$LPOHLAV[h2$person_ID == ID_01],
      age = (h2$LVEK[h2$person_ID == ID_01] - h1$Age[h1$person_ID_91 == ID_91]) %in% -1:1,
      pob = h1$BornBrno[h1$person_ID_91 == ID_91] == (h2$LMISNAR[h2$person_ID == ID_01] == 1),
      educ = h1$stvzde91_01[h1$person_ID_91 == ID_91] <= h2$LSTVZDE[h2$person_ID == ID_01]
    ) %>%
    mutate(
      join = all(sex,age,pob,educ)
    ) %>%
    ungroup() %>% 
    filter(join) %>% 
    select(ID_91,ID_01) %>% 
    distinct(ID_91, .keep_all = TRUE) %>% 
    mutate(
      adresado_01_klic = house_ID
    )
}

ivprint <- function(estm1, data = ivregX){
  require(stargazer)
  require(broom)
  require(AER)
  require(ivpack)
  
  estm1 %>% 
    lapply(summary, diagnostic = TRUE) %>% 
    lapply(magrittr::extract2, "diagnostics") %>% 
    lapply(
      function(x){
        tidy(x) %>% 
          mutate(
            statistic = round(statistic, digits = 1) %>% 
              format(digits = 1, nsmall = 1) %>% 
              str_c(.,add_stars(p.value)) %>% 
              str_trim()
          ) %>% 
          select(term = .rownames, statistic)
      }
    ) -> ivstats
  
  
  
  weak_instr <- ivstats %>% 
    lapply(function(x) filter(x, str_detect(term, "Weak")) %>% pull(statistic))
  hausman <- ivstats %>% 
    lapply(function(x) filter(x, str_detect(term, "Hausman")) %>% pull(statistic))
  sargan <- ivstats %>% 
    lapply(function(x) filter(x, str_detect(term, "Sargan")) %>% 
             mutate(statistic = ifelse(statistic == "NA","",statistic)) %>% 
             pull(statistic))
  
  
  estm1 %>% 
    lapply(function(x) cluster.robust.se(x, data$flat_ID) %>% tidy() %>% pull(std.error)) -> rse
  
  estm1 %>% 
    lapply(function(x) cluster.robust.se(x, data$flat_ID) %>% tidy() %>% pull(p.value)) -> pva
  
  stargazer(
    estm1,
    se = rse,
    p = pva,
    model.names = TRUE,
    omit = "^RCASTOBCE",
    omit.labels = "Neighborhood FE",
    type = "text",
    #dep.var.labels = "Person is unemployed",
    df = FALSE,
    add.lines = list(
      c("Weak instruments", unlist(weak_instr)),
      c("Wu-Hausman", unlist(hausman)),
      c("Sargan", unlist(sargan))
    )
  )
}

lmprint <- function(estm, data = RegData, cluster = TRUE,...){
  require(stargazer)
  require(broom)
  require(sandwich)
  require(lmtest)
  
  if(cluster){
    estm %>% 
      lapply(function(x) coeftest(x, vcov. = vcovCL(x, data$flat_ID)) %>% tidy() %>% pull(std.error)) -> rse
    
    estm %>% 
      lapply(function(x) coeftest(x, vcov. = vcovCL(x, data$flat_ID)) %>% tidy() %>% pull(p.value)) -> pva
  }else{
    estm %>% 
      lapply(function(x) coeftest(x, vcov. = vcovHC) %>% tidy() %>% pull(std.error)) -> rse
    
    estm %>% 
      lapply(function(x) coeftest(x, vcov. = vcovHC) %>% tidy() %>% pull(p.value)) -> pva
  }
  
  stargazer(
    estm,
    se = rse,
    p = pva,
    model.names = TRUE,
    omit = "^RCASTOBCE",
    omit.labels = "Neighborhood FE",
    type = "latex",
    keep.stat = c("n","adj.rsq"),
    #dep.var.labels = "Person is unemployed",
    df = FALSE,
    ...
  )
}

format_comp_table <- function(x, dig = 2){
  x %>% 
    mutate(
      mean = round(mean, digits = 3) %>% format(digits = 3, nsmall = dig) %>% str_c(.,stars) %>% str_trim(),
      sd = round(sd, digits = 3) %>% format(digits = 3, nsmall = dig, trim = TRUE) %>% str_c("(",.,")") %>% str_trim() 
    ) %>% 
    select(value,yvariable,mean,sd) %>% 
    gather(stat, coef, mean, sd) %>% 
    spread(value, coef) %>% 
    arrange(yvariable, stat)
}

get_weakinstr <- function(x){
  xs <- summary(x, diagnostics = TRUE)
  xs$diagnostic %>% 
    tidy() %>% 
    filter(str_detect(.rownames,"Weak")) %>% 
    mutate(
      statistic = round(statistic) %>% 
        #format(digits = 0, nsmall = 0) %>% 
        str_c(.,"$^{",add_stars(p.value),"}$")
    ) %>% 
    pull(statistic)
}

get_geocaliper <- function(dist_max,datasample){
  long_sd <- sd(datasample$long)
  long_mean <- mean(datasample$long)
  lat_sd <- sd(datasample$lat)
  lat_mean <- mean(datasample$lat)
  
  lat_m <- distGeo(c(long_mean,lat_mean), c(long_mean,lat_mean+lat_sd))
  long_m <- distGeo(c(long_mean,lat_mean), c(long_mean+long_sd,lat_mean))
  
  caliper_lat <- dist_max/lat_m
  caliper_long <- dist_max/long_m
  
  return(c(caliper_long, caliper_lat))
}


get_matched <- function(xmatch, caliper, geocaliper_m){
  
  geocaliper <- get_geocaliper(geocaliper_m, xmatch)
  
  suppressMessages(library(Matching))
  
  Match(
    Tr = xmatch$sample == "ATT",
    X = dplyr::select(xmatch, pscore, lat, long),
    #Weight.matrix = MM,
    caliper = c(caliper,geocaliper),
    replace = FALSE
  ) -> match
  
  data_frame(
    treatment = xmatch$adresado_01_klic[match$index.treated],
    control = xmatch$adresado_01_klic[match$index.control]
  ) %>% 
    mutate_if(is.numeric,as.integer) -> matches
  
  detach(package:Matching, unload = TRUE)
  detach(package:MASS, unload = TRUE)
  
  matches %>% 
    mutate(id = row_number()) %>% 
    gather(ctg,adresado_01_klic,-id) -> datax
  
  return(datax)
}

get_matched.g0 <- function(xmatch, caliper, geocaliper_m, tr = NULL){
  
  geocaliper <- get_geocaliper(geocaliper_m, xmatch)
  
  suppressMessages(library(Matching))
  
  if(is.null(tr)){
    Matching::Match(
      Tr = xmatch$sample == "ATT",
      X = dplyr::select(xmatch, pscore, lat, long),
      #Weight.matrix = MM,
      caliper = c(caliper,geocaliper),
      replace = TRUE
    ) -> match
  }else{
    Matching::Match(
      Tr = xmatch$sample == tr,
      X = dplyr::select(xmatch, pscore, lat, long),
      #Weight.matrix = MM,
      caliper = c(caliper,geocaliper),
      replace = TRUE
    ) -> match
  }
  
  data_frame(
    treatment = xmatch$adresado_01_klic[match$index.treated],
    control = xmatch$adresado_01_klic[match$index.control]
  ) %>% 
    mutate_if(is.numeric,as.integer) -> matches
  
  detach(package:Matching, unload = TRUE)
  detach(package:MASS, unload = TRUE)
  
  matches %>% 
    mutate(id = treatment) %>% 
    gather(ctg,adresado_01_klic, -id) %>% 
    distinct() -> datax
  
  list(
    match_obj = match,
    match = datax
  )
  
  #return(datax)
}

get_matched_NTR <- function(xmatch, caliper, geocaliper_m){
  
  geocaliper <- get_geocaliper(geocaliper_m, xmatch)
  
  suppressMessages(library(Matching))
  
  Match(
    Tr = xmatch$sample == "NTR",
    X = dplyr::select(xmatch, pscore, lat, long),
    #Weight.matrix = MM,
    caliper = c(caliper,geocaliper),
    replace = FALSE
  ) -> match
  
  data_frame(
    treatment = xmatch$adresado_01_klic[match$index.treated],
    control = xmatch$adresado_01_klic[match$index.control]
  ) %>% 
    mutate_if(is.numeric,as.integer) -> matches
  
  detach(package:Matching, unload = TRUE)
  detach(package:MASS, unload = TRUE)
  
  matches %>% 
    mutate(id = row_number()) %>% 
    gather(ctg,adresado_01_klic,-id) -> datax
  
  return(datax)
}

get_matched_RZSJ <- function(xmatch, caliper){
  
  #geocaliper <- get_geocaliper(geocaliper_m, xmatch)
  
  suppressMessages(library(Matching))
  
  Match(
    Tr = xmatch$sample == "ATT",
    X = dplyr::select(xmatch, pscore),
    #Weight.matrix = MM,
    caliper = c(caliper),
    replace = FALSE
  ) -> match
  
  data_frame(
    treatment = xmatch$adresado_01_klic[match$index.treated],
    control = xmatch$adresado_01_klic[match$index.control]
  ) %>% 
    mutate_if(is.numeric,as.integer) -> matches
  
  detach(package:Matching, unload = TRUE)
  detach(package:MASS, unload = TRUE)
  
  matches %>% 
    mutate(id = row_number()) %>% 
    gather(ctg,adresado_01_klic,-id) -> datax
  
  return(datax)
}

####

get_ownership.g0 <- function(xmatch, xdata, xmodel, clby = cluster_by){
  xdata_reg <- left_join(xmatch, xdata, by = "adresado_01_klic") #%>% 
  xdata_reg <- xdata_reg %>% drop_na()
    # group_by(person_ID) %>% 
    # mutate(
    #   WEIGHTS = 1/n() 
    # ) %>% 
    # ungroup()
  
  #emiv <- ivreg(xmodel, data = xdata_reg, weights = WEIGHTS)
  # safeIV <- safely(function(x,y) ivreg(x, data = y, weights = W))
  # safeCL <- safely(function(x,y) ivpack::cluster.robust.se(x,y) %>% tidy())
  emiv <- ivreg(xmodel, data = xdata_reg, weights = W)
  
  #emiv <- safeIV(xmodel,xdata_reg)
  
  # if(is.null(emiv$error)){
  #   emiv <- emiv$result
  # }else{
  #   out <- list(
  #     est = NA,
  #     glb = NA,
  #     coefs = NA
  #   )
  #   
  #   return(out)
  # }
  # 
  if(clby == "flat"){
    invisible(
      stop("Illegal clustering")
    #ctest <- ivpack::cluster.robust.se(emiv, xdata_reg$flat_ID) %>% tidy()
    )
  }else{
    ctest <- ivpack::cluster.robust.se(emiv, xdata_reg$adresado_01_klic) %>% tidy()
    #ctest <- (emiv, xdata_reg$adresado_01_klic)
    
  }
  
  list(
    est = emiv,
    glb = glance(emiv),
    coefs = ctest
  )
}

get_ownership.ols <- function(xmatch, xdata, xmodel, clby = cluster_by){
  xdata_reg <- left_join(xmatch, xdata, by = "adresado_01_klic") %>% drop_na()
  
  emiv <- lm(xmodel, data = xdata_reg)
  
  if(clby == "flat"){
    ctest <- coeftest(emiv, vcov. = vcovCL(emiv, xdata_reg$flat_ID)) %>% tidy() 
  }else{
    ctest <- coeftest(emiv, vcov. = vcovCL(emiv, xdata_reg$adresado_01_klic)) %>% tidy() 
  }
  
  list(
    est = emiv,
    glb = glance(emiv),
    coefs = ctest
  )
}

add_dagg <- function(x){
  y <- rep("",length(x))
  y[x >= 0.1] <- "_dagger"
  return(y)
}

# Replace RZSJ by id (matching); not used in recent version

matchingFE <- function(x){
  if(!is.Formula(x)) stop("Input object has to be of class Formula.")
  
  if("RZSJ" %in% all.vars(x)){
    if(length(x)[2] == 2){
      x <- x %>% update(.~. - factor(RZSJ) + factor(id) | . - factor(RZSJ) + factor(id))
    }else{
      x <- x %>% update(.~. - factor(RZSJ) + factor(id))
    }
  }
  
  return(x)
}


### Descriptive stats with SE

desctable <- function(x, houseCL = FALSE){
  
  if(!is.grouped_df(x)) warning("Input table is not grouped. Is it OK?")
  
  gv <- group_vars(x)
  nnmv <- select_if(x, function(x) !is.numeric(x)) %>% names()
  nnmv <- nnmv[!(nnmv %in% gv)]
  if(length(nnmv) != 0){
    str_c(nnmv, collapse = ", ") %>% str_c("Non-numeric variables: ",.," were dropped.") %>% warning()
    x <- x %>% select(-one_of(nnmv))
  }
  
  lmcl <- function(x){
    lm(value ~ 1, data = x) %>% 
      coeftest(., vcov. = vcovCL(., x$adresado_01_klic)) %>% 
      tidy()
  }
  
  
  if(houseCL){
    x %>% 
      gather(variable,value,-one_of(gv),-adresado_01_klic) %>%
      drop_na() %>% 
      group_by(variable, add = TRUE) %>%
      do(
        datar = list(.)
      ) %>% 
      ungroup() %>% 
      unnest(datar) %>% 
      group_by(sample,variable) %>% 
      do(
        reg = lmcl(.$datar[[1]])
      ) %>% 
      unnest() %>%
      select(-term) -> out
  }else{
    x %>% 
      gather(variable,value,-one_of(gv)) %>% 
      group_by(variable, add = TRUE) %>% 
      do(
        reg = lm(value ~ 1, data = .) %>% coeftest(vcov. = vcovHC) %>% tidy()
      ) %>% 
      unnest() %>%
      select(-term) -> out
  }
  
  # Number of observations
    x %>% 
      summarise(
        variable = "Observations",
        estimate = n()
      ) %>% 
      bind_rows(
        out,.
      ) -> out

  attr(out, "groups") <- gv
  
  return(out)
}

format_desctable <- function(xd, dig = 3){
  if(length(attr(xd,"groups")) != 0){
  
  xd %>% 
      rowwise() %>% 
    mutate(
      mean = estimate %>% format(digits = 1, nsmall = dig, trim = TRUE, scientific = FALSE),
      se = std.error %>% format(digits = 1, nsmall = dig, trim = TRUE, scientific = FALSE) %>% str_c("(",.,")")
    ) %>% 
      ungroup() %>% 
    select(one_of(attr(xd,"groups")),variable,mean,se) %>% 
    gather(stat,value,-variable,-one_of(attr(xd,"groups"))) %>%
    unite(sample,one_of(attr(xd,"groups"))) %>% 
    spread(sample,value) %>% 
    arrange(variable,stat) -> out
    
  }else{
    
    xd %>% 
      rowwise() %>% 
      mutate(
        mean = estimate %>% format(digits = 1, nsmall = dig, trim = TRUE, scientific = FALSE),
        se = std.error %>% format(digits = 1, nsmall = dig, trim = TRUE, scientific = FALSE) %>% str_c("(",.,")")
      ) %>% 
      ungroup() %>% 
      select(variable,mean,se) %>% 
      gather(stat,value,-variable) %>%
      arrange(variable,stat) -> out
    
  }
  
  out %>% 
    filter(!(variable == "Observations" & stat == "se")) %>% 
    mutate(
      stat = ifelse(variable == "Observations","n",stat)
    )
}


# Print

add_stars <- function(x, latex = FALSE, ds = FALSE){
  y <- character(length(x))
  
  if(latex){
    if(ds){
      y[x<0.1] <- "$^{*}$"
      y[x<0.05] <- "$^{**}$"
      y[x<0.01] <- "$^{***}$"
    }else{
      y[x<0.1] <- "^{*}"
      y[x<0.05] <- "^{**}"
      y[x<0.01] <- "^{***}"
    }
  }else{
    y[x<0.1] <- "*"
    y[x<0.05] <- "**"
    y[x<0.01] <- "***"
  }
  
  return(y)
}

print_latex <- function(x, path){
  write.table(x, 
              file = path, 
              quote = FALSE, 
              sep = "\t&\t", 
              eol = "\\\\\n", 
              col.names = FALSE, 
              row.names = FALSE,
              na = ""
  )
}

get_reg_column <- function(x, est_d = 3, se_d = 3, dsr = FALSE){
  out <- x %>% 
    filter(!str_detect(term,"RZSJ")) %>% 
    filter(!str_detect(term,"LVEK")) %>% 
    group_by(term) %>% 
    transmute(
      est = estimate %>% #round(digits = est_d) %>% 
        format(digits = 1, nsmall = est_d, trim = TRUE, scientific = FALSE) %>% 
        str_c(.,add_stars(p.value, latex = TRUE, ds = dsr)),
      se = std.error %>% #round(digits = se_d) %>% 
        format(digits = 1, nsmall = se_d, trim = TRUE, scientific = FALSE) %>% 
        str_c("(",.,")")
    ) %>% 
    ungroup() %>% 
    gather(stat,value,-term)
  
  tribble(
    ~term, ~stat, ~value,
    "LVEKfe", "fe", ifelse(any(str_detect(x$term,"LVEK")),"\\c{Yes}","\\c{--}"),
    "RZSJfe", "fe", ifelse(any(str_detect(x$term,"RZSJ")),"\\c{Yes}","\\c{--}")
  ) %>% 
    bind_rows(
      out,.
    )
}

get_tables <- function(estm,estm_ctest, ffile = NULL, space = TRUE, ivtest = FALSE, 
                       houses = FALSE, labels = NULL, reorder = NULL, intercept = "auto",
                       dollar = FALSE, dig = 3){
  
  
  if(houses){
    flev <- c(
      "EDUCmiddle",
      "EDUChigh",
      "LVEK",
      "MALE",
      "BornBrno",
      "byty[10,20)",
      "byty[20,40)",
      "(Intercept)"
    )
    
    flab <- c(
      "Secondary education (share)",
      "Tertiary education (share)",
      "Age (mean)",
      "Male (share)",
      "Born in Brno (share)",
      "No. of apartments $\\in [10,20)$",
      "No. of apartments $\\geq 20$",
      "Constant"
    )
  }else{
    flev <- c(
      "ownershipowner",
      "sample_dummyTRUE",
      "sample_dummy",
      "I(sample != \"NTR\")TRUE",
      "EDUCmiddle",
      "EDUChigh",
      "EXP",
      "I((EXP^2)/100)",
      "LPOHLAVmale",
      "BornBrno",
      "htypeIncomplete_family",
      "htypeOne_person_household",
      "(Intercept)",
      "LVEKfe",
      "RZSJfe"
    )
    
    flab <- c(
      "Homeowner (=1)",
      "Living in privatized house (=1)",
      "Living in privatized house (=1)",
      "Living in city-owned house (=1)",
      "Secondary education (=1)",
      "Tertiary education (=1)",
      "Potential experience",
      "Potential experience$^2/100$",
      "Male (=1)",
      "Born in Brno (=1)",
      "Incomplete family (=1)",
      "One-person family (=1)",
      "Constant",
      "Age dummies",
      "Neighborhood dummies"
    )
  }
  
  if(!is.null(labels)){
    flev <- labels$flev
    flab <- labels$flab
  }
  
  if(!is.null(reorder)){
    estm <- estm[reorder]
    estm_ctest <- estm_ctest[reorder]
  }
  
  
  if(intercept == "auto"){
    intercept <- estm %>% 
      lapply(function(x) names(x$coefficients)) %>% 
      lapply(function(x){
        c(
          any(str_detect(x,"RZSJ")),
          any(str_detect(x,"LVEK"))
        ) %>% any()
        }) %>% 
      unlist()
  }
  
  #### FE
  #estm %>% lapply(function(x) names(x$coefficients) %>% str_detect("LVEK") %>% any()) %>% unlist()
  #estm %>% lapply(function(x) names(x$coefficients) %>% str_detect("RZSJ") %>% any()) %>% unlist()
  
  r2 <- estm %>% lapply(class) %>% unlist() %>% unique()
  r2 <- !("glm" %in% r2)
  
  estm_ctest %>%
    lapply(tidy) %>% 
    lapply(function(x) x %>% mutate(term = ifelse(str_detect(term,"sample_dummy"),"sample_dummy",term))) %>% 
    lapply(get_reg_column, dsr = dollar) %>% 
    reduce(full_join, by = c("term","stat")) -> aux
    
  if(all(intercept)){
    aux <- aux %>% filter(!str_detect(term,"Intercept"))
  }else if(any(intercept)){
    intercept <- c(FALSE,FALSE,intercept)
    
    in_rows <- aux %>% 
      filter(!str_detect(term,"Intercept"))
    
    i_rows <- aux %>% 
      filter(str_detect(term,"Intercept"))
    
    i_rows[,intercept] <- NA
    
    aux <- bind_rows(
      in_rows,
      i_rows
    )
  }
  
  aux %>% 
    mutate(
      term = ifelse(str_detect(term,"ownership"),"ownershipowner",term)
    ) %>% 
    mutate(
      term = term %>% factor(
        levels = flev,
        labels = flab
      )
    ) %>% 
    arrange(term,stat) %>% 
    mutate_if(is.factor,as.character) %>% 
    mutate_all(as.character) %>% 
    mutate(
      term = ifelse(stat == "se",NA,term)
    ) %>% 
    select(-stat) -> b1
  
  if(space){
    b1 <- b1 %>% 
      mutate(
        rown = row_number(),
        term = ifelse(
          rown == 1,
          term,
          ifelse(
            !is.na(term) & rown != 1,
            str_c("\\addlinespace\n ",term),
            NA
          )
        )
      ) %>% 
      select(-rown)
  }
  
  if(ivtest){
  estm %>% 
    lapply(summary, diagnostic = TRUE) %>% 
    lapply(magrittr::extract2, "diagnostics") %>% 
    lapply(as_tibble) %>% 
    lapply(function(x){
      x[1,3] <- x[1,3]/100
      return(x)
      }) %>% 
    lapply(function(x) transmute(x, out = #round(statistic, digits = 3) %>% 
                                   format(statistic, digits = 1, nsmall = dig, trim = TRUE, scientific = FALSE) %>% 
                                   str_c(.,add_stars(`p-value`, ds = dollar)) %>% str_trim())) -> ivdiag
  
  #weak_test <- ivdiag %>% lapply(function(x) magrittr::extract2(x,"out")[1]) %>% unlist
  #hausman_test <- ivdiag %>% lapply(function(x) magrittr::extract2(x,"out")[2]) %>% unlist
  #sargan_test <- ivdiag %>% lapply(function(x) magrittr::extract2(x,"out")[3]) %>% unlist
  
  data_frame(
    rname = letters[1:length(estm)],
    r2 = estm %>% lapply(glance) %>% bind_rows() %>% 
      rowwise() %>% 
      mutate(
        r.squared = r.squared %>% format(nsmall=dig, digits = 1, scientific = FALSE, trim = TRUE)
      ) %>% pull(r.squared),
    #round(digits = est_d) %>% format(nsmall = est_d, trim = TRUE),
    obs = estm %>% lapply(nobs) %>% unlist() %>% as.integer() %>% format(trim = TRUE, big.mark = ","),# %>% str_c("\\c{",.,"}"),
    #weak = weak_test,
    #hausman = hausman_test
  ) %>% 
    gather(stat,value,-rname) %>% 
    mutate(
      value = str_c("\\c{",value,"}")
    ) %>% 
    spread(rname,value) %>% 
    mutate(
      stat = stat %>% factor(levels = c("weak","hausman","obs","r2"))
    ) %>% 
    arrange(stat) %>% 
    select(-stat) %>% 
    bind_cols(
      data_frame(
        #term = c("Weak instruments /100","Wu-Hausman","Observations", "R^2")
        term = c("Observations", "$R^2$")
      ),
      .
    ) -> b2
  }else{
  
    if(r2){
  data_frame(
    rname = letters[1:length(estm)],
    r2 = estm %>% lapply(glance) %>% bind_rows() %>% 
      rowwise() %>% 
      mutate(
        r.squared = r.squared %>% format(nsmall=dig, digits = 1, scientific = FALSE, trim = TRUE)
      ) %>% pull(r.squared),
      #round(digits = est_d) %>% format(nsmall = est_d, trim = TRUE),
    obs = estm %>% lapply(nobs) %>% unlist() %>% as.integer() %>% format(trim = TRUE, big.mark = ",") #%>% str_c("\\c{",.,"}")
  ) %>% 
    gather(stat,value,-rname) %>% 
    mutate(
      value = str_c("\\c{",value,"}")
    ) %>% 
    spread(rname,value) %>% 
    arrange(stat) %>% 
    select(-stat) %>% 
    bind_cols(
      data_frame(
        term = c("Observations", "$R^2$")
      ),
      .
    ) -> b2
    }else{
      data_frame(
        rname = letters[1:length(estm)],
        # r2 = estm %>% lapply(glance) %>% bind_rows() %>% 
        #   rowwise() %>% 
        #   mutate(
        #     adj.r.squared = adj.r.squared %>% format(nsmall=3, digits = 1, scientific = FALSE, trim = TRUE)
        #   ) %>% pull(adj.r.squared),
        #round(digits = est_d) %>% format(nsmall = est_d, trim = TRUE),
        obs = estm %>% lapply(nobs) %>% unlist() %>% as.integer() %>% format(trim = TRUE, big.mark = ",") %>% str_c("\\c{",.,"}")
      ) %>% 
        gather(stat,value,-rname) %>% 
        mutate(
          value = str_c("\\c{",value,"}")
        ) %>% 
        spread(rname,value) %>% 
        arrange(stat) %>% 
        select(-stat) %>% 
        bind_cols(
          data_frame(
            term = c("Observations")
          ),
          .
        ) -> b2
    }
  }
  
  if(!is.null(ffile)){
    print_latex(b1,str_c(ffile,"_b1.tex"))
    print_latex(b2,str_c(ffile,"_b2.tex"))
  }else{
    names(b1) <- names(b2)
    bind_rows(b1,b2) %>% 
      mutate_all(str_remove_all,"\\\\c") %>% 
      mutate_all(str_remove_all,"\\{") %>% 
      mutate_all(str_remove_all,"\\}") %>% 
      mutate_all(str_remove_all,"\\^")
  }
}

matching_geo_restrictions <- function(sampleH,geoMAX,preload = NULL){
  
  if(is.null(preload)){
    gzsj <- st_read("GIS/AdministrativniCleneni_v13.gdb/","ZakladniSidelniJednotkyBody", quiet = TRUE)
  }else{
    gzsj <- preload
  }
  
  gzsj <- gzsj %>% filter(NAZ_OBEC == "Brno")
  DIST <- st_distance(gzsj) %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    magrittr::set_colnames(gzsj$KOD_ZSJ) %>% 
    mutate_all(as.double) %>% 
    mutate(
      RZSJ_1 = gzsj$KOD_ZSJ
    ) %>% 
    gather(RZSJ_2,distM,-RZSJ_1)
  
  
  rzsj_list <- sampleH %>% 
    mutate(rn = row_number()) %>% 
    select(rn,RZSJ)
  
  crossing(
    RN_1 = rzsj_list$rn,
    RN_2 = rzsj_list$rn
  ) %>% 
    left_join(
      rzsj_list, by = c("RN_1"="rn")
    ) %>% 
    rename(
      RZSJ_1 = RZSJ
    ) %>% 
    left_join(
      rzsj_list, by = c("RN_2"="rn")
    ) %>% 
    rename(
      RZSJ_2 = RZSJ
    ) %>% 
    left_join(DIST, by = c("RZSJ_1", "RZSJ_2")) %>%
    filter(RN_1 != RN_2) %>% 
    filter(distM > geoMAX) %>%
    mutate(notallowed = -1) %>% 
    select(RN_1,RN_2,notallowed) %>% 
    as.matrix()
}

get_matched.geo <- function(xmatch, caliper, maxdist, tr = NULL, geo = NULL){
  
  if(is.null(tr)){
    Matching::Match(
      Tr = xmatch$sample == "ATT",
      X = dplyr::select(xmatch, pscore),
      #Weight.matrix = MM,
      caliper = c(caliper),
      ties = TRUE,
      restrict = matching_geo_restrictions(xmatch,maxdist,preload = geo),
      replace = TRUE
    ) -> match
  }else{
    Matching::Match(
      Tr = xmatch$sample == tr,
      X = dplyr::select(xmatch, pscore),
      #Weight.matrix = MM,
      caliper = c(caliper),
      ties = TRUE,
      restrict = matching_geo_restrictions(xmatch,maxdist,preload = geo),
      replace = TRUE
    ) -> match
  }
  
  data_frame(
    treatment = xmatch$adresado_01_klic[match$index.treated],
    control = xmatch$adresado_01_klic[match$index.control]
  ) %>% 
    mutate_if(is.numeric,as.integer) -> matches
  
  #detach(package:Matching, unload = TRUE)
  #detach(package:MASS, unload = TRUE)
  
  matches %>% 
    mutate(id = treatment) %>% 
    group_by(treatment) %>% 
    add_tally() %>% 
    ungroup() %>% 
    mutate(
      W = 1/n
    ) %>% 
    select(-n) %>% 
    gather(ctg,adresado_01_klic, -id, -W) %>% 
    distinct() -> datax
  
  list(
    match_obj = match,
    match = datax
  )
}

is_equal <- function(x,y, tolerance = 1e-4){
  if(!is.numeric(x) | !is.numeric(y)) stop("Inputs must be numeric.")
  if(length(y) != 1 & length(x)!=length(y)) stop("y must be of length 1 or equal length of x")
  
  abs(x-y) < tolerance
}
