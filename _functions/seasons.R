seasons <- function(year, location) {
  
  # Requires these packages
  require(lutz)
  require(insol)
  require(sf)
  
  # Example inputs
  # x <- st_read(file.path(AOI_dir, "0_raw_inputs", "base_layers", "TRIM", "aoi.gpkg")) %>% 
  #   st_centroid()
  # year <- 2013
  
  ## Goal: Calculate the season of a given location/date input
  ## problem: insol package rounds values. Need to load the package, but then
  ## edit the daylength function to remove the rounding problem
  
  # Required functions
  # From insol package
  daylength2 <- function (lat, long, jd, tmz) {
    if (nargs() < 4 ) {cat("USAGE: daylength(latitude, longitude, jd, timezone) \n values in degrees, julian days, hours \n"); return()}
    EqTime = eqtime(jd)
    delta = declination(jd)
    tanlatdel = -tan(radians(lat)) * tan(radians(delta))
    tanlatdel[tanlatdel>1]=1
    omega = acos(tanlatdel)
    daylen = (2*omega)/(2*pi/24)
    stndmeridian = tmz*15
    deltaLatTime = long-stndmeridian
    deltaLatTime = deltaLatTime * 24/360 
    sunrise = 12*(1-omega/pi)-deltaLatTime-EqTime/60 
    sunset = 12*(1+omega/pi)-deltaLatTime-EqTime/60
    sunrise[omega==0] = NA
    sunset[omega==0] = NA
    return(cbind(sunrise,sunset,daylen))
  }
  
  # From fmdates package
  # julian_day_to_gregorian <- function(julian_day) {
  #   # Jean-Meeus Astronomical Algorithms, Chapter 7, pg 63
  #   z <- trunc(julian_day + 0.5)
  #   f <- julian_day + 0.5 - z
  #   alpha <- trunc((z - 1867216.25) / 36524.25)
  #   a <- rep(NA, NROW(julian_day))
  #   a[z < 2299161] <- z
  #   a[!(z < 2299161)] <- z + 1 + alpha - trunc(alpha / 4)
  #   b <- a + 1524
  #   c <- trunc((b - 122.1) / 365.25)
  #   d <- trunc(365.25 * c)
  #   e <- trunc((b - d) / 30.6001)
  #   dom <- b - d - trunc(30.6001 * e) + f
  #   f_dom <- dom - trunc(dom)
  #   m1 <- e - 1
  #   m2 <- rep(0, NROW(e))
  #   m2[e == 14 | e == 15] <-  (e - 13)[e == 14 | e == 15]
  #   m <- m1 * (e < 14) + m2 * (e >= 14)
  #   m1 <- c - 4716
  #   m2 <- rep(0, NROW(m))
  #   m2[m == 1 | m == 2] <- (c - 4715)[m == 1 | m == 2]
  #   y <- m1 * (m > 2) + m2 * (m <= 2)
  #   # reproduce same output as default values of lubridate::ymd()
  #   ISOdate(y, m, trunc(dom), hour = 0, tz = 'UTC') + f_dom * 24 * 60 * 60
  # }
  # 
  # # From fmdates package
  # delta_t <- function(dates) {
  #   assertthat::assert_that(lubridate::is.POSIXt(dates))
  #   dates <- lubridate::with_tz(dates, "UTC")
  #   # See: http://maia.usno.navy.mil/ser7/tai-utc.dat
  #   # Ten leap seconds before start of data contained in .leap.seconds
  #   # http://www.usno.navy.mil/USNO/earth-orientation/eo-info/general/date-time-def/date-and-time-definitions
  #   n_leap_seconds <- Map(function (x)
  #     lubridate::with_tz(.leap.seconds, "UTC") <= x, dates)
  #   n_leap_seconds <- vapply(n_leap_seconds, sum, integer(1))
  #   # Observed:
  #   # http://maia.usno.navy.mil/ser7/deltat.data
  #   # Quick comparison to observed suggests this code is correct to within a
  #   # second or two.
  #   32.184 + (10 + n_leap_seconds)
  # }
  # 
  # # From fmdates package
  # jde_to_gregorian <- function(jde, tz = "UTC", want_dt = FALSE) {
  #   res <- julian_day_to_gregorian(jde)
  #   if (!want_dt) {
  #     # Meeus pg 71 has Delta T = TD - UT.
  #     res <- lubridate::floor_date(res - lubridate::seconds(delta_t(res)),
  #                                  unit = "minute")
  #     
  #   }
  #   if (tz != "UTC") {
  #     return (lubridate::with_tz(res, tz = tz))
  #   } else {
  #     return (res)
  #   }
  # }
  # 
  # # From fmdates package
  # rads <- function(degs) {
  #   degs * pi / 180
  # }
  # 
  # # From fmdates package
  # equinox <- function(years, season = 'mar', tz = "UTC", want_dt = FALSE) {
  #   # Check inputs
  #   assertthat::assert_that(season %in% c('mar', 'sep'),
  #                           all(years >= 1000), all(years <= 3000))
  #   
  #   # Algorithm from Jean-Meeus Astronomical Algorithms, Chapter 26
  #   # Table 26.B coefficients
  #   coef_mar <- c(2451623.80984, 365242.37404,  0.05169, -0.00411, -0.00057)
  #   coef_jun <- c(2451716.56767, 365241.62603,  0.00325, 0.00888, -0.00030)
  #   coef_sep <- c(2451810.21715, 365242.01767, -0.11575,  0.00337,  0.00078)
  #   
  #   # Table 26.C constants.
  #   A <- c(485, 203, 199, 182, 156, 136, 77, 74, 70, 58, 52, 50, 45, 44, 29, 18,
  #          17, 16, 14, 12, 12, 12, 9, 8)
  #   B <- c(324.96, 337.23, 342.08, 27.85, 73.14, 171.52, 222.54, 296.72, 243.58,
  #          119.81, 297.17, 21.02, 247.54, 325.15, 60.93, 155.12, 288.79, 198.04,
  #          199.76, 95.39, 287.11, 320.81, 227.73, 15.45)
  #   C <- c(1934.136, 32964.467, 20.186, 445267.112, 45036.886, 22518.443,
  #          65928.934, 3034.906, 9037.513, 33718.147, 150.678, 2281.226,
  #          29929.562, 31555.956, 4443.417, 67555.328, 4562.452, 62894.029,
  #          31436.921, 14577.848, 31931.756, 34777.259, 1222.114, 16859.074)
  #   
  #   # Mean time
  #   y <- (years - 2000) / 1000
  #   M <- cbind(rep(1, NROW(y)), y, y ^ 2, y ^ 3, y ^ 4)
  #   if (identical(season, 'mar')) {
  #     jde0 <- as.vector(M %*% coef_mar)
  #   } else if (identical(season, "jun")) {
  #     jde0 <- as.vector(M %*% coef_jun)
  #   } else {
  #     jde0 <- as.vector(M %*% coef_sep)
  #   }
  #   # Correction
  #   tt <- (jde0 - 2451545) / 36525
  #   w <- 35999.373 * tt - 2.47
  #   delta_lambda <- 1 + 0.0334 * cos(rads(w)) + 0.0007 * cos(rads(2 * w))
  #   s <- vector("numeric", length(tt))
  #   for (i in seq_along(s)) {
  #     s[i] <- Reduce(sum, Map(function (a, b, c, t) a * cos(rads(b + c * t)),
  #                             A, B, C, tt[i]))
  #   }
  #   # julian datetime in dynamical time (i.e. not UTC)
  #   jde <- jde0 + 0.00001 * s / delta_lambda
  #   jde_to_gregorian(jde, tz, want_dt)
  # }
  
  ## With functions set, calculate the dates that seasons start
  
  x_tz <- tz_list() %>% 
    dplyr::filter(tz_name == tz_lookup(location, method = "accurate"))
  
  dates <- seq(ISOdate(year, 1, 1, tz = x_tz$tz_name[1]), 
               ISOdate(year, 12, 31, tz = x_tz$tz_name[1]), by = 'day')
  j_dates <- JD(dates)
  
  # Winter overlaps years, so get next years winter dates
  dates_next <- seq(ISOdate(year + 1, 1, 1, tz = x_tz$tz_name[1]), 
               ISOdate(year + 1, 12, 31, tz = x_tz$tz_name[1]), by = 'day')
  j_dates_next <- JD(dates_next)
  
  coords <- st_transform(location, 4326L) %>% st_coordinates()
  
  dl <- as.data.frame(daylength2(lat = coords[2], long = coords[1], jd = j_dates, tmz = x_tz$utc_offset_h[1])) %>% 
    mutate(date = dates)
  dl_next <- as.data.frame(daylength2(lat = coords[2], long = coords[1], jd = j_dates_next, tmz = x_tz$utc_offset_h[1])) %>% 
    mutate(date = dates_next)
  dl_1 <- slice_head(dl, n = ceiling(nrow(dl) / 2))
  dl_2 <- slice_tail(dl, n = floor(nrow(dl) / 2))
  dl_1_next <- slice_head(dl_next, n = ceiling(nrow(dl_next) / 2))
  
  spring <- slice_min(dl_1, abs(daylen - 12)) %>% dplyr::pull(date)
  summer <- slice_max(dl, daylen) %>% dplyr::pull(date)
  fall <- slice_min(dl_2, abs(daylen - 12)) %>% dplyr::pull(date)
  winter <- slice_min(dl, daylen) %>% dplyr::pull(date)
  spring_next <- slice_min(dl_1_next, abs(daylen - 12)) %>% dplyr::pull(date)
  
  # spring <- equinox(year, season = "mar", tz = x_tz$tz_name[1])
  # fall <- equinox(year, season = "sep", tz = x_tz$tz_name[1])
  
  seasons <- data.frame(
    season = c("spring", "summer", "fall", "winter"),
    start = c(spring, summer, fall, winter),
    end = c(summer, fall, winter, spring_next)) %>% 
    dplyr::mutate(across(c(start, end), as.Date)) %>% 
    dplyr::mutate(end = end - 1)
  
}
