#calculate stuff for jb
library(chillR)
library(tidyverse)


flist <- list.files('future_weather/', full.names = TRUE)

future_weather <- purrr::map(flist, function(file_path) {
  
  # Extract filename
  filename <- basename(file_path) # just the file name, no directory
  
  # Split by '.' and extract parts
  name_parts <- str_split(filename, '\\.')[[1]]
  
  # Check if name_parts is long enough
  if (length(name_parts) < 5) {
    stop(paste("Filename", filename, "does not have enough parts."))
  }
  
  # Read and process the file
  read.csv(file_path) %>%
    select(-X) %>%
    mutate(
      ssp = name_parts[3],
      gcm = name_parts[4],
      scenario_year = name_parts[5],
      yday = lubridate::yday(DATE),
      season = ifelse(yday >= 175, Year + 1, Year)
    )
})

future_weather <- do.call('rbind', future_weather) 
future_weather$id = paste(future_weather$ssp, future_weather$gcm, future_weather$scenario_year, sep = '--')


#weed out incomplete seasons
drop_list <- future_weather %>% 
  group_by(id, season) %>% 
  summarise(n = n()) %>% 
  filter(n < 365)

future_weather <- future_weather %>% 
  filter(!(paste(id, season) %in% paste(drop_list$id, drop_list$season)))

hist_weather <- read.csv('weather_2020_koeln-bonn.csv') %>% 
  mutate(scenario_year = 2020, 
         ssp = 'historical',
         gcm = 'historical',
         yday = lubridate::yday(DATE),
         season = ifelse(yday >= 175,
                         yes = Year  +1,
                         no = Year),
         id = paste(ssp, gcm, scenario_year, sep = '--'))
drop_list <- hist_weather %>% 
  group_by(id, season) %>% 
  summarise(n = n()) %>% 
  filter(n < 365)
hist_weather <- hist_weather %>% 
  filter(!(paste(id, season) %in% paste(drop_list$id, drop_list$season)))

#combine hist and future weather
weather_combined <- future_weather %>% 
  rbind(hist_weather)


weather_combined <- weather_combined %>%
  group_by(ssp, season) %>%
  mutate(id_season = cur_group_id() - 1) %>%
  ungroup()


#split by each season
test <- split(x = weather_combined, f = weather_combined$id_seaon) 


#get ssp names
ssp <- names(test) %>% 
  str_split('--') %>% 
  purrr::map_chr(1)

#split by ssp
future_weather <- split(x = future_weather, f = future_weather$ssp)

example_season <- future_weather %>% 
  filter(season == 2001, ssp == 'ssp126',
         gcm == 'ACCESS-CM2')



example_season <- example_season %>% 
  mutate(yday_plot = ifelse(yday >= 175, yes = yday - yday_subtract, no = yday)) 


write.csv(example_season, file = 'example_season.csv', row.names = FALSE)






# #threshold what is considered to be rain
# rain_cutoff <- 1
# alpha_soil <- 0.2
# alpha_wet_reduce <- 0.8
# alpha_wet <- alpha_soil * alpha_wet_reduce
# 
# #calculate soil temperature etc
# example_season <- example_season %>% 
#   mutate(soil_wet = check_soil_wet(Prec),
#          Tmean = (Tmin + Tmax) / 2,
#          alpha = ifelse(soil_wet, yes = alpha_wet, no = alpha_soil),
#          T_soil = get_Tsoil(Tmean, alpha = alpha)) 
# 
# #plot evolution of temperatures
# example_season %>% 
# ggplot(aes(x = yday_plot)) +
#   geom_line(aes(y = Tmax, col = 'Tmax')) +
#   geom_line(aes(y = Tmin, col = 'Tmin')) +
#   geom_line(aes(y = T_soil, col = 'T_soil')) +
#   geom_col(aes(y = Prec, fill = 'Prec')) +
#   scale_fill_manual(values = 'grey50') +
#   theme_bw()
# ggsave('example_temperature.jpeg', height = 10, width = 15,
#        unit = 'cm', device = 'jpeg')








#------------------------#
#period since harvest
#------------------------#

weather_harvest <- weather_sub[i_havest:nrow(weather_sub),]




