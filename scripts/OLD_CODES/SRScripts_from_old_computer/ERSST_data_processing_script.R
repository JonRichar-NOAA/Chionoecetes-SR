#####################################################################################################################################
#################################################     ERSST SST     #################################################################
#################################################import data set and define vars#####################################################

# Extract Bering Sea SST data from netCDF file of global monthly temperatures
# on 2x2 degree grid estimated by interpolation of available data by NOAA's Climate Data Center  


library(ncdf4)
library(chron)
library(maps)
library(mapdata)
library(lmtest)          ####for durbin-watson AR test
library(nlme)            ##### for gls models with autoregression
library(lattice)

##########################################################################################################################################
########################################## Attach and define data ########################################################################

# Open netCDF file:
ersst <- nc_open("data/sst.mnmean.V4.nc")



########################################### Check the dataset:############################################################

ersst    # basic info about dimensions and variables
lon <- ncvar_get(ersst,varid="lon")
lon    # view longitudes (degrees East) - every 2 degrees

lat <- ncvar_get(ersst,varid="lon")     
lat
# view latitudes - every 2 degrees

# view dates (middle of month):
get.var.ncdf(ersst, "time")  # Days since January 1, 1800 
d <- dates(get.var.ncdf(ersst, "time"), origin=c(1,15,1800)) #converts dates to mid-month (15th) for every month since 1954 to 2009
d 


#d[1504]
#d[1516]
#d[1528]
#d[1540]
#d[1552]
#d[1564]
#d[1576]
#d[1588]
#d[1600]
#d[1612]
#d[1624]
#d[1636]
#d[1648]
#d[1660]
#d[1672]
#d[1684]
#d[1696]
#d[1708]
#d[1720]
#d[1732]
#d[1744]
#d[1756]
#d[1768]
#d[1780]
#d[1792]
#d[1804]
#d[1816]

SST <- ncvar_get(ersst,varid="sst") 
sst<-na.omit(SST)

################################## 1978 ########################################
April.78<-sst[95:101,15:18,1492]
May.78<-sst[95:101,15:18,1493]
June.78<-sst[95:101,15:18,1494]
July.78<-sst[95:101,15:18,1495]
April.May.78<-sst[95:101,15:18,1492:1493]
May.June.78<-sst[95:101,15:18,1493:1494]
June.July.78<-sst[95:101,15:18,1494:1495]
April.July.78<-sst[95:101,15:18,1492:1495]
May.July.78<-sst[95:101,15:18,1493:1495]

mean.April.78<-mean(April.78, na.rm=T)
mean.May.78<-mean(May.78,na.rm=T)
mean.June.78<-mean(June.78,na.rm=T)
mean.July.78<-mean(July.78,na.rm=T)
mean.April.May.78<-mean(April.May.78,na.rm=T)
mean.May.June.78<-mean(May.June.78,na.rm=T)
mean.June.July.78<-mean(June.July.78,na.rm=T)
mean.April.July.78<-mean(April.July.78,na.rm=T)
mean.May.July.78<-mean(May.July.78,na.rm=T)

################################## 1979 ########################################
April.79<-sst[95:101,15:18,1504]
May.79<-sst[95:101,15:18,1505]
June.79<-sst[95:101,15:18,1506]
July.79<-sst[95:101,15:18,1507]
April.May.79<-sst[95:101,15:18,1504:1505]
May.June.79<-sst[95:101,15:18,1505:1506]
June.July.79<-sst[95:101,15:18,1506:1507]
April.July.79<-sst[95:101,15:18,1504:1507]
May.July.79<-sst[95:101,15:18,1505:1507]


mean.April.79<-mean(April.79, na.rm=T)
mean.May.79<-mean(May.79,na.rm=T)
mean.June.79<-mean(June.79,na.rm=T)
mean.July.79<-mean(July.79,na.rm=T)
mean.April.May.79<-mean(April.May.79,na.rm=T)
mean.May.June.79<-mean(May.June.79,na.rm=T)
mean.June.July.79<-mean(June.July.79,na.rm=T)
mean.April.July.79<-mean(April.July.79,na.rm=T)
mean.May.July.79<-mean(May.July.79,na.rm=T)

################################## 1980 ########################################
April.80<-sst[95:101,15:18,1516]
May.80<-sst[95:101,15:18,1517]
June.80<-sst[95:101,15:18,1518]
July.80<-sst[95:101,15:18,1519]
April.May.80<-sst[95:101,15:18,1516:1517]
May.June.80<-sst[95:101,15:18,1517:1518]
June.July.80<-sst[95:101,15:18,1518:1519]
April.July.80<-sst[95:101,15:18,1516:1519]
May.July.80<-sst[95:101,15:18,1517:1519]

mean.April.80<-mean(April.80, na.rm=T)
mean.May.80<-mean(May.80,na.rm=T)
mean.June.80<-mean(June.80,na.rm=T)
mean.July.80<-mean(July.80,na.rm=T)
mean.April.May.80<-mean(April.May.80,na.rm=T)
mean.May.June.80<-mean(May.June.80,na.rm=T)
mean.June.July.80<-mean(June.July.80,na.rm=T)
mean.April.July.80<-mean(April.July.80,na.rm=T)
mean.May.July.80<-mean(May.July.80,na.rm=T)

################################## 1981 ########################################
April.81<-sst[95:101,15:18,1528]
May.81<-sst[95:101,15:18,1529]
June.81<-sst[95:101,15:18,1530]
July.81<-sst[95:101,15:18,1531]
April.May.81<-sst[95:101,15:18,1528:1529]
May.June.81<-sst[95:101,15:18,1529:1530]
June.July.81<-sst[95:101,15:18,1530:1531]
April.July.81<-sst[95:101,15:18,1528:1531]
May.July.81<-sst[95:101,15:18,1529:1531]


mean.April.81<-mean(April.81, na.rm=T)
mean.May.81<-mean(May.81,na.rm=T)
mean.June.81<-mean(June.81,na.rm=T)
mean.July.81<-mean(July.81,na.rm=T)
mean.April.May.81<-mean(April.May.81,na.rm=T)
mean.May.June.81<-mean(May.June.81,na.rm=T)
mean.June.July.81<-mean(June.July.81,na.rm=T)
mean.April.July.81<-mean(April.July.81,na.rm=T)
mean.May.July.81<-mean(May.July.81,na.rm=T)

################################## 1982 ########################################
April.82<-sst[95:101,15:18,1540]
May.82<-sst[95:101,15:18,1541]
June.82<-sst[95:101,15:18,1542]
July.82<-sst[95:101,15:18,1543]
April.May.82<-sst[95:101,15:18,1540:1541]
May.June.82<-sst[95:101,15:18,1541:1542]
June.July.82<-sst[95:101,15:18,1542:1543]
April.July.82<-sst[95:101,15:18,1540:1543]
May.July.82<-sst[95:101,15:18,1541:1543]


mean.April.82<-mean(April.82, na.rm=T)
mean.May.82<-mean(May.82,na.rm=T)
mean.June.82<-mean(June.82,na.rm=T)
mean.July.82<-mean(July.82,na.rm=T)
mean.April.May.82<-mean(April.May.82,na.rm=T)
mean.May.June.82<-mean(May.June.82,na.rm=T)
mean.June.July.82<-mean(June.July.82,na.rm=T)
mean.April.July.82<-mean(April.July.82,na.rm=T)
mean.May.July.82<-mean(May.July.82,na.rm=T)

################################## 1983 ########################################
April.83<-sst[95:101,15:18,1552]
May.83<-sst[95:101,15:18,1553]
June.83<-sst[95:101,15:18,1554]
July.83<-sst[95:101,15:18,1555]
April.May.83<-sst[95:101,15:18,1552:1553]
May.June.83<-sst[95:101,15:18,1553:1554]
June.July.83<-sst[95:101,15:18,1554:1555]
April.July.83<-sst[95:101,15:18,1552:1555]
May.July.83<-sst[95:101,15:18,1553:1555]


mean.April.83<-mean(April.83, na.rm=T)
mean.May.83<-mean(May.83,na.rm=T)
mean.June.83<-mean(June.83,na.rm=T)
mean.July.83<-mean(July.83,na.rm=T)
mean.April.May.83<-mean(April.May.83,na.rm=T)
mean.May.June.83<-mean(May.June.83,na.rm=T)
mean.June.July.83<-mean(June.July.83,na.rm=T)
mean.April.July.83<-mean(April.July.83,na.rm=T)
mean.May.July.83<-mean(May.July.83,na.rm=T)

################################## 1984 ########################################
April.84<-sst[95:101,15:18,1564]
May.84<-sst[95:101,15:18,1565]
June.84<-sst[95:101,15:18,1566]
July.84<-sst[95:101,15:18,1567]
April.May.84<-sst[95:101,15:18,1564:1565]
May.June.84<-sst[95:101,15:18,1565:1566]
June.July.84<-sst[95:101,15:18,1566:1567]
April.July.84<-sst[95:101,15:18,1564:1567]
May.July.84<-sst[95:101,15:18,1565:1567]


mean.April.84<-mean(April.84, na.rm=T)
mean.May.84<-mean(May.84,na.rm=T)
mean.June.84<-mean(June.84,na.rm=T)
mean.July.84<-mean(July.84,na.rm=T)
mean.April.May.84<-mean(April.May.84,na.rm=T)
mean.May.June.84<-mean(May.June.84,na.rm=T)
mean.June.July.84<-mean(June.July.84,na.rm=T)
mean.April.July.84<-mean(April.July.84,na.rm=T)
mean.May.July.84<-mean(May.July.84,na.rm=T)


################################## 1985 ########################################
April.85<-sst[95:101,15:18,1576]
May.85<-sst[95:101,15:18,1577]
June.85<-sst[95:101,15:18,1578]
July.85<-sst[95:101,15:18,1579]
April.May.85<-sst[95:101,15:18,1576:1577]
May.June.85<-sst[95:101,15:18,1577:1578]
June.July.85<-sst[95:101,15:18,1578:1579]
April.July.85<-sst[95:101,15:18,1576:1579]
May.July.85<-sst[95:101,15:18,1577:1579]


mean.April.85<-mean(April.85, na.rm=T)
mean.May.85<-mean(May.85,na.rm=T)
mean.June.85<-mean(June.85,na.rm=T)
mean.July.85<-mean(July.85,na.rm=T)
mean.April.May.85<-mean(April.May.85,na.rm=T)
mean.May.June.85<-mean(May.June.85,na.rm=T)
mean.June.July.85<-mean(June.July.85,na.rm=T)
mean.April.July.85<-mean(April.July.85,na.rm=T)
mean.May.July.85<-mean(May.July.85,na.rm=T)

################################## 1986 ########################################
April.86<-sst[95:101,15:18,1588]
May.86<-sst[95:101,15:18,1589]
June.86<-sst[95:101,15:18,1590]
July.86<-sst[95:101,15:18,1591]
April.May.86<-sst[95:101,15:18,1588:1589]
May.June.86<-sst[95:101,15:18,1589:1590]
June.July.86<-sst[95:101,15:18,1590:1591]
April.July.86<-sst[95:101,15:18,1588:1591]
May.July.86<-sst[95:101,15:18,1589:1591]


mean.April.86<-mean(April.86, na.rm=T)
mean.May.86<-mean(May.86,na.rm=T)
mean.June.86<-mean(June.86,na.rm=T)
mean.July.86<-mean(July.86,na.rm=T)
mean.April.May.86<-mean(April.May.86,na.rm=T)
mean.May.June.86<-mean(May.June.86,na.rm=T)
mean.June.July.86<-mean(June.July.86,na.rm=T)
mean.April.July.86<-mean(April.July.86,na.rm=T)
mean.May.July.86<-mean(May.July.86,na.rm=T)


################################## 1987 ########################################
April.87<-sst[95:101,15:18,1600]
May.87<-sst[95:101,15:18,1601]
June.87<-sst[95:101,15:18,1602]
July.87<-sst[95:101,15:18,1603]
April.May.87<-sst[95:101,15:18,1600:1601]
May.June.87<-sst[95:101,15:18,1601:1602]
June.July.87<-sst[95:101,15:18,1602:1603]
April.July.87<-sst[95:101,15:18,1600:1603]
May.July.87<-sst[95:101,15:18,1601:1603]


mean.April.87<-mean(April.87, na.rm=T)
mean.May.87<-mean(May.87,na.rm=T)
mean.June.87<-mean(June.87,na.rm=T)
mean.July.87<-mean(July.87,na.rm=T)
mean.April.May.87<-mean(April.May.87,na.rm=T)
mean.May.June.87<-mean(May.June.87,na.rm=T)
mean.June.July.87<-mean(June.July.87,na.rm=T)
mean.April.July.87<-mean(April.July.87,na.rm=T)
mean.May.July.87<-mean(May.July.87,na.rm=T)


################################## 1988 ########################################
April.88<-sst[95:101,15:18,1612]
May.88<-sst[95:101,15:18,1613]
June.88<-sst[95:101,15:18,1614]
July.88<-sst[95:101,15:18,1615]
April.May.88<-sst[95:101,15:18,1612:1613]
May.June.88<-sst[95:101,15:18,1613:1614]
June.July.88<-sst[95:101,15:18,1614:1615]
April.July.88<-sst[95:101,15:18,1612:1615]
May.July.88<-sst[95:101,15:18,1613:1615]


mean.April.88<-mean(April.88, na.rm=T)
mean.May.88<-mean(May.88,na.rm=T)
mean.June.88<-mean(June.88,na.rm=T)
mean.July.88<-mean(July.88,na.rm=T)
mean.April.May.88<-mean(April.May.88,na.rm=T)
mean.May.June.88<-mean(May.June.88,na.rm=T)
mean.June.July.88<-mean(June.July.88,na.rm=T)
mean.April.July.88<-mean(April.July.88,na.rm=T)
mean.May.July.88<-mean(May.July.88,na.rm=T)


################################## 1989 ########################################
April.89<-sst[95:101,15:18,1624]
May.89<-sst[95:101,15:18,1625]
June.89<-sst[95:101,15:18,1626]
July.89<-sst[95:101,15:18,1627]
April.May.89<-sst[95:101,15:18,1624:1625]
May.June.89<-sst[95:101,15:18,1625:1626]
June.July.89<-sst[95:101,15:18,1626:1627]
April.July.89<-sst[95:101,15:18,1624:1627]
May.July.89<-sst[95:101,15:18,1625:1627]


mean.April.89<-mean(April.89, na.rm=T)
mean.May.89<-mean(May.89,na.rm=T)
mean.June.89<-mean(June.89,na.rm=T)
mean.July.89<-mean(July.89,na.rm=T)
mean.April.May.89<-mean(April.May.89,na.rm=T)
mean.May.June.89<-mean(May.June.89,na.rm=T)
mean.June.July.89<-mean(June.July.89,na.rm=T)
mean.April.July.89<-mean(April.July.89,na.rm=T)
mean.May.July.89<-mean(May.July.89,na.rm=T)


################################## 1990 ########################################
April.90<-sst[95:101,15:18,1636]
May.90<-sst[95:101,15:18,1637]
June.90<-sst[95:101,15:18,1638]
July.90<-sst[95:101,15:18,1639]
April.May.90<-sst[95:101,15:18,1636:1637]
May.June.90<-sst[95:101,15:18,1637:1638]
June.July.90<-sst[95:101,15:18,1638:1639]
April.July.90<-sst[95:101,15:18,1636:1639]
May.July.90<-sst[95:101,15:18,1637:1639]


mean.April.90<-mean(April.90, na.rm=T)
mean.May.90<-mean(May.90,na.rm=T)
mean.June.90<-mean(June.90,na.rm=T)
mean.July.90<-mean(July.90,na.rm=T)
mean.April.May.90<-mean(April.May.90,na.rm=T)
mean.May.June.90<-mean(May.June.90,na.rm=T)
mean.June.July.90<-mean(June.July.90,na.rm=T)
mean.April.July.90<-mean(April.July.90,na.rm=T)
mean.May.July.90<-mean(May.July.90,na.rm=T)



################################## 1991 ########################################
April.91<-sst[95:101,15:18,1648]
May.91<-sst[95:101,15:18,1649]
June.91<-sst[95:101,15:18,1650]
July.91<-sst[95:101,15:18,1651]
April.May.91<-sst[95:101,15:18,1648:1649]
May.June.91<-sst[95:101,15:18,1649:1650]
June.July.91<-sst[95:101,15:18,1650:1651]
April.July.91<-sst[95:101,15:18,1648:1651]
May.July.91<-sst[95:101,15:18,1649:1651]



mean.April.91<-mean(April.91, na.rm=T)
mean.May.91<-mean(May.91,na.rm=T)
mean.June.91<-mean(June.91,na.rm=T)
mean.July.91<-mean(July.91,na.rm=T)
mean.April.May.91<-mean(April.May.91,na.rm=T)
mean.May.June.91<-mean(May.June.91,na.rm=T)
mean.June.July.91<-mean(June.July.91,na.rm=T)
mean.April.July.91<-mean(April.July.91,na.rm=T)
mean.May.July.91<-mean(May.July.91,na.rm=T)


################################## 1992 ########################################
April.92<-sst[95:101,15:18,1660]
May.92<-sst[95:101,15:18,1661]
June.92<-sst[95:101,15:18,1662]
July.92<-sst[95:101,15:18,1663]
April.May.92<-sst[95:101,15:18,1660:1661]
May.June.92<-sst[95:101,15:18,1661:1662]
June.July.92<-sst[95:101,15:18,1662:1663]
April.July.92<-sst[95:101,15:18,1660:1663]
May.July.92<-sst[95:101,15:18,1661:1663]


mean.April.92<-mean(April.92, na.rm=T)
mean.May.92<-mean(May.92,na.rm=T)
mean.June.92<-mean(June.92,na.rm=T)
mean.July.92<-mean(July.92,na.rm=T)
mean.April.May.92<-mean(April.May.92,na.rm=T)
mean.May.June.92<-mean(May.June.92,na.rm=T)
mean.June.July.92<-mean(June.July.92,na.rm=T)
mean.April.July.92<-mean(April.July.92,na.rm=T)
mean.May.July.92<-mean(May.July.92,na.rm=T)


################################## 1993 ########################################
April.93<-sst[95:101,15:18,1672]
May.93<-sst[95:101,15:18,1673]
June.93<-sst[95:101,15:18,1674]
July.93<-sst[95:101,15:18,1675]
April.May.93<-sst[95:101,15:18,1672:1673]
May.June.93<-sst[95:101,15:18,1673:1674]
June.July.93<-sst[95:101,15:18,1674:1675]
April.July.93<-sst[95:101,15:18,1672:1675]
May.July.93<-sst[95:101,15:18,1673:1675]


mean.April.93<-mean(April.93, na.rm=T)
mean.May.93<-mean(May.93,na.rm=T)
mean.June.93<-mean(June.93,na.rm=T)
mean.July.93<-mean(July.93,na.rm=T)
mean.April.May.93<-mean(April.May.93,na.rm=T)
mean.May.June.93<-mean(May.June.93,na.rm=T)
mean.June.July.93<-mean(June.July.93,na.rm=T)
mean.April.July.93<-mean(April.July.93,na.rm=T)
mean.May.July.93<-mean(May.July.93,na.rm=T)

################################## 1994 ########################################
April.94<-sst[95:101,15:18,1684]
May.94<-sst[95:101,15:18,1685]
June.94<-sst[95:101,15:18,1686]
July.94<-sst[95:101,15:18,1687]
April.May.94<-sst[95:101,15:18,1684:1685]
May.June.94<-sst[95:101,15:18,1685:1686]
June.July.94<-sst[95:101,15:18,1686:1687]
April.July.94<-sst[95:101,15:18,1684:1687]
May.July.94<-sst[95:101,15:18,1685:1687]

mean.April.94<-mean(April.94, na.rm=T)
mean.May.94<-mean(May.94,na.rm=T)
mean.June.94<-mean(June.94,na.rm=T)
mean.July.94<-mean(July.94,na.rm=T)
mean.April.May.94<-mean(April.May.94,na.rm=T)
mean.May.June.94<-mean(May.June.94,na.rm=T)
mean.June.July.94<-mean(June.July.94,na.rm=T)
mean.April.July.94<-mean(April.July.94,na.rm=T)
mean.May.July.94<-mean(May.July.94,na.rm=T)


################################## 1995 ########################################
April.95<-sst[95:101,15:18,1696]
May.95<-sst[95:101,15:18,1697]
June.95<-sst[95:101,15:18,1698]
July.95<-sst[95:101,15:18,1699]
April.May.95<-sst[95:101,15:18,1696:1697]
May.June.95<-sst[95:101,15:18,1697:1698]
June.July.95<-sst[95:101,15:18,1698:1699]
April.July.95<-sst[95:101,15:18,1696:1699]
May.July.95<-sst[95:101,15:18,1697:1699]


mean.April.95<-mean(April.95, na.rm=T)
mean.May.95<-mean(May.95,na.rm=T)
mean.June.95<-mean(June.95,na.rm=T)
mean.July.95<-mean(July.95,na.rm=T)
mean.April.May.95<-mean(April.May.95,na.rm=T)
mean.May.June.95<-mean(May.June.95,na.rm=T)
mean.June.July.95<-mean(June.July.95,na.rm=T)
mean.April.July.95<-mean(April.July.95,na.rm=T)
mean.May.July.95<-mean(May.July.95,na.rm=T)

################################## 1996 ########################################
April.96<-sst[95:101,15:18,1708]
May.96<-sst[95:101,15:18,1709]
June.96<-sst[95:101,15:18,1710]
July.96<-sst[95:101,15:18,1711]
April.May.96<-sst[95:101,15:18,1708:1709]
May.June.96<-sst[95:101,15:18,1709:1710]
June.July.96<-sst[95:101,15:18,1710:1711]
April.July.96<-sst[95:101,15:18,1708:1711]
May.July.96<-sst[95:101,15:18,1709:1711]

mean.April.96<-mean(April.96, na.rm=T)
mean.May.96<-mean(May.96,na.rm=T)
mean.June.96<-mean(June.96,na.rm=T)
mean.July.96<-mean(July.96,na.rm=T)
mean.April.May.96<-mean(April.May.96,na.rm=T)
mean.May.June.96<-mean(May.June.96,na.rm=T)
mean.June.July.96<-mean(June.July.96,na.rm=T)
mean.April.July.96<-mean(April.July.96,na.rm=T)
mean.May.July.96<-mean(May.July.96,na.rm=T)


################################## 1997 ########################################
April.97<-sst[95:101,15:18,1720]
May.97<-sst[95:101,15:18,1721]
June.97<-sst[95:101,15:18,1722]
July.97<-sst[95:101,15:18,1723]
April.May.97<-sst[95:101,15:18,1720:1721]
May.June.97<-sst[95:101,15:18,1721:1722]
June.July.97<-sst[95:101,15:18,1722:1723]
April.July.97<-sst[95:101,15:18,1720:1723]
May.July.97<-sst[95:101,15:18,1721:1723]


mean.April.97<-mean(April.97, na.rm=T)
mean.May.97<-mean(May.97,na.rm=T)
mean.June.97<-mean(June.97,na.rm=T)
mean.July.97<-mean(July.97,na.rm=T)
mean.April.May.97<-mean(April.May.97,na.rm=T)
mean.May.June.97<-mean(May.June.97,na.rm=T)
mean.June.July.97<-mean(June.July.97,na.rm=T)
mean.April.July.97<-mean(April.July.97,na.rm=T)
mean.May.July.97<-mean(May.July.97,na.rm=T)

################################## 1998 ########################################
April.98<-sst[95:101,15:18,1732]
May.98<-sst[95:101,15:18,1733]
June.98<-sst[95:101,15:18,1734]
July.98<-sst[95:101,15:18,1735]
April.May.98<-sst[95:101,15:18,1732:1733]
May.June.98<-sst[95:101,15:18,1733:1734]
June.July.98<-sst[95:101,15:18,1734:1735]
April.July.98<-sst[95:101,15:18,1732:1735]
May.July.98<-sst[95:101,15:18,1733:1735]

mean.April.98<-mean(April.98, na.rm=T)
mean.May.98<-mean(May.98,na.rm=T)
mean.June.98<-mean(June.98,na.rm=T)
mean.July.98<-mean(July.98,na.rm=T)
mean.April.May.98<-mean(April.May.98,na.rm=T)
mean.May.June.98<-mean(May.June.98,na.rm=T)
mean.June.July.98<-mean(June.July.98,na.rm=T)
mean.April.July.98<-mean(April.July.98,na.rm=T)
mean.May.July.98<-mean(May.July.98,na.rm=T)



################################## 1999 ########################################
April.99<-sst[95:101,15:18,1744]
May.99<-sst[95:101,15:18,1745]
June.99<-sst[95:101,15:18,1746]
July.99<-sst[95:101,15:18,1747]
April.May.99<-sst[95:101,15:18,1744:1745]
May.June.99<-sst[95:101,15:18,1745:1746]
June.July.99<-sst[95:101,15:18,1746:1747]
April.July.99<-sst[95:101,15:18,1744:1747]
May.July.99<-sst[95:101,15:18,1745:1747]


mean.April.99<-mean(April.99, na.rm=T)
mean.May.99<-mean(May.99,na.rm=T)
mean.June.99<-mean(June.99,na.rm=T)
mean.July.99<-mean(July.99,na.rm=T)
mean.April.May.99<-mean(April.May.99,na.rm=T)
mean.May.June.99<-mean(May.June.99,na.rm=T)
mean.June.July.99<-mean(June.July.99,na.rm=T)
mean.April.July.99<-mean(April.July.99,na.rm=T)
mean.May.July.99<-mean(May.July.99,na.rm=T)

################################## 2000 ########################################
April.00<-sst[95:101,15:18,1756]
May.00<-sst[95:101,15:18,1757]
June.00<-sst[95:101,15:18,1758]
July.00<-sst[95:101,15:18,1759]
April.May.00<-sst[95:101,15:18,1756:1757]
May.June.00<-sst[95:101,15:18,1757:1758]
June.July.00<-sst[95:101,15:18,1758:1759]
April.July.00<-sst[95:101,15:18,1756:1759]
May.July.00<-sst[95:101,15:18,1757:1759]


mean.April.00<-mean(April.00, na.rm=T)
mean.May.00<-mean(May.00,na.rm=T)
mean.June.00<-mean(June.00,na.rm=T)
mean.July.00<-mean(July.00,na.rm=T)
mean.April.May.00<-mean(April.May.00,na.rm=T)
mean.May.June.00<-mean(May.June.00,na.rm=T)
mean.June.July.00<-mean(June.July.00,na.rm=T)
mean.April.July.00<-mean(April.July.00,na.rm=T)
mean.May.July.00<-mean(May.July.00,na.rm=T)


################################## 2001 ########################################
April.01<-sst[95:101,15:18,1768]
May.01<-sst[95:101,15:18,1769]
June.01<-sst[95:101,15:18,1770]
July.01<-sst[95:101,15:18,1771]
April.May.01<-sst[95:101,15:18,1768:1769]
May.June.01<-sst[95:101,15:18,1769:1770]
June.July.01<-sst[95:101,15:18,1770:1771]
April.July.01<-sst[95:101,15:18,1768:1771]
May.July.01<-sst[95:101,15:18,1769:1771]


mean.April.01<-mean(April.01, na.rm=T)
mean.May.01<-mean(May.01,na.rm=T)
mean.June.01<-mean(June.01,na.rm=T)
mean.July.01<-mean(July.01,na.rm=T)
mean.April.May.01<-mean(April.May.01,na.rm=T)
mean.May.June.01<-mean(May.June.01,na.rm=T)
mean.June.July.01<-mean(June.July.01,na.rm=T)
mean.April.July.01<-mean(April.July.01,na.rm=T)
mean.May.July.01<-mean(May.July.01,na.rm=T)


################################## 2002 ########################################
April.02<-sst[95:101,15:18,1780]
May.02<-sst[95:101,15:18,1781]
June.02<-sst[95:101,15:18,1782]
July.02<-sst[95:101,15:18,1783]
April.May.02<-sst[95:101,15:18,1780:1781]
May.June.02<-sst[95:101,15:18,1781:1782]
June.July.02<-sst[95:101,15:18,1782:1783]
April.July.02<-sst[95:101,15:18,1780:1783]
May.July.02<-sst[95:101,15:18,1781:1783]


mean.April.02<-mean(April.02, na.rm=T)
mean.May.02<-mean(May.02,na.rm=T)
mean.June.02<-mean(June.02,na.rm=T)
mean.July.02<-mean(July.02,na.rm=T)
mean.April.May.02<-mean(April.May.02,na.rm=T)
mean.May.June.02<-mean(May.June.02,na.rm=T)
mean.June.July.02<-mean(June.July.02,na.rm=T)
mean.April.July.02<-mean(April.July.02,na.rm=T)
mean.May.July.02<-mean(May.July.02,na.rm=T)


################################## 2003 ########################################
April.03<-sst[95:101,15:18,1792]
May.03<-sst[95:101,15:18,1793]
June.03<-sst[95:101,15:18,1794]
July.03<-sst[95:101,15:18,1795]
April.May.03<-sst[95:101,15:18,1792:1793]
May.June.03<-sst[95:101,15:18,1793:1794]
June.July.03<-sst[95:101,15:18,1794:1795]
April.July.03<-sst[95:101,15:18,1792:1795]
May.July.03<-sst[95:101,15:18,1793:1795]


mean.April.03<-mean(April.03, na.rm=T)
mean.May.03<-mean(May.03,na.rm=T)
mean.June.03<-mean(June.03,na.rm=T)
mean.July.03<-mean(July.03,na.rm=T)
mean.April.May.03<-mean(April.May.03,na.rm=T)
mean.May.June.03<-mean(May.June.03,na.rm=T)
mean.June.July.03<-mean(June.July.03,na.rm=T)
mean.April.July.03<-mean(April.July.03,na.rm=T)
mean.May.July.03<-mean(May.July.03,na.rm=T)


################################## 2004 ########################################
April.04<-sst[95:101,15:18,1804]
May.04<-sst[95:101,15:18,1805]
June.04<-sst[95:101,15:18,1806]
July.04<-sst[95:101,15:18,1807]
April.May.04<-sst[95:101,15:18,1804:1805]
May.June.04<-sst[95:101,15:18,1805:1806]
June.July.04<-sst[95:101,15:18,1806:1807]
April.July.04<-sst[95:101,15:18,1804:1807]
May.July.04<-sst[95:101,15:18,1805:1807]


mean.April.04<-mean(April.04, na.rm=T)
mean.May.04<-mean(May.04,na.rm=T)
mean.June.04<-mean(June.04,na.rm=T)
mean.July.04<-mean(July.04,na.rm=T)
mean.April.May.04<-mean(April.May.04,na.rm=T)
mean.May.June.04<-mean(May.June.04,na.rm=T)
mean.June.July.04<-mean(June.July.04,na.rm=T)
mean.April.July.04<-mean(April.July.04,na.rm=T)
mean.May.July.04<-mean(May.July.04,na.rm=T)



################################## 2005 ########################################
April.05<-sst[95:101,15:18,1816]
May.05<-sst[95:101,15:18,1817]
June.05<-sst[95:101,15:18,1818]
July.05<-sst[95:101,15:18,1819]
April.May.05<-sst[95:101,15:18,1816:1817]
May.June.05<-sst[95:101,15:18,1817:1818]
June.July.05<-sst[95:101,15:18,1818:1819]
April.July.05<-sst[95:101,15:18,1816:1819]
May.July.05<-sst[95:101,15:18,1817:1819]


mean.April.05<-mean(April.05, na.rm=T)
mean.May.05<-mean(May.05,na.rm=T)
mean.June.05<-mean(June.05,na.rm=T)
mean.July.05<-mean(July.05,na.rm=T)
mean.April.May.05<-mean(April.May.05,na.rm=T)
mean.May.June.05<-mean(May.June.05,na.rm=T)
mean.June.July.05<-mean(June.July.05,na.rm=T)
mean.April.July.05<-mean(April.July.05,na.rm=T)
mean.May.July.05<-mean(May.July.05,na.rm=T)



################################## 2006 ########################################
April.06<-sst[95:101,15:18,1828]
May.06<-sst[95:101,15:18,1829]
June.06<-sst[95:101,15:18,1830]
July.06<-sst[95:101,15:18,1831]
April.May.06<-sst[95:101,15:18,1828:1829]
May.June.06<-sst[95:101,15:18,1829:1830]
June.July.06<-sst[95:101,15:18,1830:1831]
April.July.06<-sst[95:101,15:18,1828:1831]
May.July.06<-sst[95:101,15:18,1829:1831]


mean.April.06<-mean(April.06, na.rm=T)
mean.May.06<-mean(May.06,na.rm=T)
mean.June.06<-mean(June.06,na.rm=T)
mean.July.06<-mean(July.06,na.rm=T)
mean.April.May.06<-mean(April.May.06,na.rm=T)
mean.May.June.06<-mean(May.June.06,na.rm=T)
mean.June.July.06<-mean(June.July.06,na.rm=T)
mean.April.July.06<-mean(April.July.06,na.rm=T)
mean.May.July.06<-mean(May.July.06,na.rm=T)

################################## 2007 ########################################
April.07<-sst[95:101,15:18,1840]
May.07<-sst[95:101,15:18,1841]
June.07<-sst[95:101,15:18,1842]
July.07<-sst[95:101,15:18,1843]
April.May.07<-sst[95:101,15:18,1840:1841]
May.June.07<-sst[95:101,15:18,1841:1842]
June.July.07<-sst[95:101,15:18,1842:1843]
April.July.07<-sst[95:101,15:18,1840:1843]
May.July.07<-sst[95:101,15:18,1841:1843]


mean.April.07<-mean(April.07, na.rm=T)
mean.May.07<-mean(May.07,na.rm=T)
mean.June.07<-mean(June.07,na.rm=T)
mean.July.07<-mean(July.07,na.rm=T)
mean.April.May.07<-mean(April.May.07,na.rm=T)
mean.May.June.07<-mean(May.June.07,na.rm=T)
mean.June.July.07<-mean(June.July.07,na.rm=T)
mean.April.July.07<-mean(April.July.07,na.rm=T)
mean.May.July.07<-mean(May.July.07,na.rm=T)

################################## 2008 ########################################
April.08<-sst[95:101,15:18,1840]
May.08<-sst[95:101,15:18,1841]
June.08<-sst[95:101,15:18,1842]
July.08<-sst[95:101,15:18,1843]
April.May.08<-sst[95:101,15:18,1840:1841]
May.June.08<-sst[95:101,15:18,1841:1842]
June.July.08<-sst[95:101,15:18,1842:1843]
April.July.08<-sst[95:101,15:18,1840:1843]
May.July.08<-sst[95:101,15:18,1841:1843]


mean.April.08<-mean(April.08, na.rm=T)
mean.May.08<-mean(May.08,na.rm=T)
mean.June.08<-mean(June.08,na.rm=T)
mean.July.08<-mean(July.08,na.rm=T)
mean.April.May.08<-mean(April.May.08,na.rm=T)
mean.May.June.08<-mean(May.June.08,na.rm=T)
mean.June.July.08<-mean(June.July.08,na.rm=T)
mean.April.July.08<-mean(April.July.08,na.rm=T)
mean.May.July.08<-mean(May.July.08,na.rm=T)

###########9####################### 2009 ########################################
April.09<-sst[95:101,15:18,1852]
May.09<-sst[95:101,15:18,1853]
June.09<-sst[95:101,15:18,1854]
July.09<-sst[95:101,15:18,1855]
April.May.09<-sst[95:101,15:18,1852:1853]
May.June.09<-sst[95:101,15:18,1853:1854]
June.July.09<-sst[95:101,15:18,1854:1855]
April.July.09<-sst[95:101,15:18,1852:1855]
May.July.09<-sst[95:101,15:18,1853:1855]


mean.April.09<-mean(April.09, na.rm=T)
mean.May.09<-mean(May.09,na.rm=T)
mean.June.09<-mean(June.09,na.rm=T)
mean.July.09<-mean(July.09,na.rm=T)
mean.April.May.09<-mean(April.May.09,na.rm=T)
mean.May.June.09<-mean(May.June.09,na.rm=T)
mean.June.July.09<-mean(June.July.09,na.rm=T)
mean.April.July.09<-mean(April.July.09,na.rm=T)
mean.May.July.09<-mean(May.July.09,na.rm=T)

###########9####################### 2010 ########################################
April.10<-sst[95:101,15:18,1864]
May.10<-sst[95:101,15:18,1865]
June.10<-sst[95:101,15:18,1866]
July.10<-sst[95:101,15:18,1867]
April.May.10<-sst[95:101,15:18,1864:1865]
May.June.10<-sst[95:101,15:18,1865:1866]
June.July.10<-sst[95:101,15:18,1866:1867]
April.July.10<-sst[95:101,15:18,1864:1867]
May.July.10<-sst[95:101,15:18,1865:1867]


mean.April.10<-mean(April.10, na.rm=T)
mean.May.10<-mean(May.10,na.rm=T)
mean.June.10<-mean(June.10,na.rm=T)
mean.July.10<-mean(July.10,na.rm=T)
mean.April.May.10<-mean(April.May.10,na.rm=T)
mean.May.June.10<-mean(May.June.10,na.rm=T)
mean.June.July.10<-mean(June.July.10,na.rm=T)
mean.April.July.10<-mean(April.July.10,na.rm=T)
mean.May.July.10<-mean(May.July.10,na.rm=T)

###########9####################### 2011 ########################################
April.11<-sst[95:101,15:18,1876]
May.11<-sst[95:101,15:18,1877]
June.11<-sst[95:101,15:18,1878]
July.11<-sst[95:101,15:18,1879]
April.May.11<-sst[95:101,15:18,1876:1877]
May.June.11<-sst[95:101,15:18,1877:1878]
June.July.11<-sst[95:101,15:18,1878:1879]
April.July.11<-sst[95:101,15:18,1876:1879]
May.July.11<-sst[95:101,15:18,1877:1879]


mean.April.11<-mean(April.11, na.rm=T)
mean.May.11<-mean(May.11,na.rm=T)
mean.June.11<-mean(June.11,na.rm=T)
mean.July.11<-mean(July.11,na.rm=T)
mean.April.May.11<-mean(April.May.11,na.rm=T)
mean.May.June.11<-mean(May.June.11,na.rm=T)
mean.June.July.11<-mean(June.July.11,na.rm=T)
mean.April.July.11<-mean(April.July.11,na.rm=T)
mean.May.July.11<-mean(May.July.11,na.rm=T)

###########9####################### 2012 ########################################
April.12<-sst[95:101,15:18,1876]
May.12<-sst[95:101,15:18,1877]
June.12<-sst[95:101,15:18,1878]
July.12<-sst[95:101,15:18,1879]
April.May.12<-sst[95:101,15:18,1876:1877]
May.June.12<-sst[95:101,15:18,1877:1878]
June.July.12<-sst[95:101,15:18,1878:1879]
April.July.12<-sst[95:101,15:18,1876:1879]
May.July.12<-sst[95:101,15:18,1877:1879]


mean.April.12<-mean(April.12, na.rm=T)
mean.May.12<-mean(May.12,na.rm=T)
mean.June.12<-mean(June.12,na.rm=T)
mean.July.12<-mean(July.12,na.rm=T)
mean.April.May.12<-mean(April.May.12,na.rm=T)
mean.May.June.12<-mean(May.June.12,na.rm=T)
mean.June.July.12<-mean(June.July.12,na.rm=T)
mean.April.July.12<-mean(April.July.12,na.rm=T)
mean.May.July.12<-mean(May.July.12,na.rm=T)

###########9####################### 2013 ########################################
April.13<-sst[95:101,15:18,1888]
May.13<-sst[95:101,15:18,1889]
June.13<-sst[95:101,15:18,1890]
July.13<-sst[95:101,15:18,1891]
April.May.13<-sst[95:101,15:18,1888:1889]
May.June.13<-sst[95:101,15:18,1889:1890]
June.July.13<-sst[95:101,15:18,1890:1891]
April.July.13<-sst[95:101,15:18,1888:1891]
May.July.13<-sst[95:101,15:18,1889:1891]


mean.April.13<-mean(April.13, na.rm=T)
mean.May.13<-mean(May.13,na.rm=T)
mean.June.13<-mean(June.13,na.rm=T)
mean.July.13<-mean(July.13,na.rm=T)
mean.April.May.13<-mean(April.May.13,na.rm=T)
mean.May.June.13<-mean(May.June.13,na.rm=T)
mean.June.July.13<-mean(June.July.13,na.rm=T)
mean.April.July.13<-mean(April.July.13,na.rm=T)
mean.May.July.13<-mean(May.July.13,na.rm=T)


###########9####################### 2014 ########################################
April.14<-sst[95:101,15:18,1900]
May.14<-sst[95:101,15:18,1901]
June.14<-sst[95:101,15:18,1902]
July.14<-sst[95:101,15:18,1903]
April.May.14<-sst[95:101,15:18,1900:1901]
May.June.14<-sst[95:101,15:18,1901:1902]
June.July.14<-sst[95:101,15:18,1902:1903]
April.July.14<-sst[95:101,15:18,1900:1903]
May.July.14<-sst[95:101,15:18,1901:1903]


mean.April.14<-mean(April.14, na.rm=T)
mean.May.14<-mean(May.14,na.rm=T)
mean.June.14<-mean(June.14,na.rm=T)
mean.July.14<-mean(July.14,na.rm=T)
mean.April.May.14<-mean(April.May.14,na.rm=T)
mean.May.June.14<-mean(May.June.14,na.rm=T)
mean.June.July.14<-mean(June.July.14,na.rm=T)
mean.April.July.14<-mean(April.July.14,na.rm=T)
mean.May.July.14<-mean(May.July.14,na.rm=T)

###########9####################### 2015 ########################################
April.15<-sst[95:101,15:18,1912]
May.15<-sst[95:101,15:18,1913]
June.15<-sst[95:101,15:18,1914]
July.15<-sst[95:101,15:18,1915]
April.May.15<-sst[95:101,15:18,1912:1913]
May.June.15<-sst[95:101,15:18,1913:1914]
June.July.15<-sst[95:101,15:18,1914:1915]
April.July.15<-sst[95:101,15:18,1912:1915]
May.July.15<-sst[95:101,15:18,1913:1915]


mean.April.15<-mean(April.15, na.rm=T)
mean.May.15<-mean(May.15,na.rm=T)
mean.June.15<-mean(June.15,na.rm=T)
mean.July.15<-mean(July.15,na.rm=T)
mean.April.May.15<-mean(April.May.15,na.rm=T)
mean.May.June.15<-mean(May.June.15,na.rm=T)
mean.June.July.15<-mean(June.July.15,na.rm=T)
mean.April.July.15<-mean(April.July.15,na.rm=T)
mean.May.July.15<-mean(May.July.15,na.rm=T)

###########9####################### 2016 ########################################
April.16<-sst[95:101,15:18,1924]
May.16<-sst[95:101,15:18,1925]
June.16<-sst[95:101,15:18,1926]
July.16<-sst[95:101,15:18,1927]
April.May.16<-sst[95:101,15:18,1924:1925]
May.June.16<-sst[95:101,15:18,1925:1926]
June.July.16<-sst[95:101,15:18,1926:1927]
April.July.16<-sst[95:101,15:18,1924:1927]
May.July.16<-sst[95:101,15:18,1925:1927]


mean.April.16<-mean(April.16, na.rm=T)
mean.May.16<-mean(May.16,na.rm=T)
mean.June.16<-mean(June.16,na.rm=T)
mean.July.16<-mean(July.16,na.rm=T)
mean.April.May.16<-mean(April.May.16,na.rm=T)
mean.May.June.16<-mean(May.June.16,na.rm=T)
mean.June.July.16<-mean(June.July.16,na.rm=T)
mean.April.July.16<-mean(April.July.16,na.rm=T)
mean.May.July.16<-mean(May.July.16,na.rm=T)
#####################################################################################################################################
################################################### Combine into vectors for analysis ###############################################

April.sst<-c(mean.April.78,mean.April.79,mean.April.80,mean.April.81,mean.April.82,mean.April.83,mean.April.84,mean.April.85,mean.April.86,
             mean.April.87,mean.April.88,mean.April.89,mean.April.90,mean.April.91,mean.April.92,mean.April.93,mean.April.94,mean.April.95,mean.April.96,
             mean.April.97,mean.April.98,mean.April.99,mean.April.00,mean.April.01,mean.April.02,mean.April.03,mean.April.04,mean.April.05,
             mean.April.06,mean.April.07,mean.April.08,mean.April.09,mean.April.10,mean.April.11,mean.April.12,mean.April.13,mean.April.14,
             mean.April.15,mean.April.16)

May.sst<-c(mean.May.78,mean.May.79,mean.May.80,mean.May.81,mean.May.82,mean.May.83,mean.May.84,mean.May.85,mean.May.86,mean.May.87,mean.May.88,
           mean.May.89,mean.May.90,mean.May.91,mean.May.92,mean.May.93,mean.May.94,mean.May.95,mean.May.96,mean.May.97,mean.May.98,mean.May.99,mean.May.00,
           mean.May.01,mean.May.02,mean.May.03,mean.May.04,mean.May.05,mean.May.06,mean.May.07,mean.May.08,mean.May.09,mean.May.10,
           mean.May.11,mean.May.12,mean.May.13,mean.May.14,mean.May.15,mean.May.16)

June.sst<-c(mean.June.78,mean.June.79,mean.June.80,mean.June.81,mean.June.82,mean.June.83,mean.June.84,mean.June.85,mean.June.86,mean.June.87,
            mean.June.88,mean.June.89,mean.June.90,mean.June.91,mean.June.92,mean.June.93,mean.June.94,mean.June.95,mean.June.96,mean.June.97,mean.June.98,
            mean.June.99,mean.June.00,mean.June.01,mean.June.02,mean.June.03,mean.June.04,mean.June.05,mean.June.06,mean.June.07,mean.June.08,
            mean.June.09,mean.June.10,mean.June.11,mean.June.12,mean.June.13,mean.June.14,mean.June.15,mean.June.16)

July.sst<-c(mean.July.78,mean.July.79,mean.July.80,mean.July.81,mean.July.82,mean.July.83,mean.July.84,mean.July.85,mean.July.86,mean.July.87,
            mean.July.88,mean.July.89,mean.July.90,mean.July.91,mean.July.92,mean.July.93,mean.July.94,mean.July.95,mean.July.96,mean.July.97,mean.July.98,
            mean.July.99,mean.July.00,mean.July.01,mean.July.02,mean.July.03,mean.July.04,mean.July.05,mean.July.06,mean.July.07,mean.July.08,
            mean.July.09,mean.July.10,mean.July.11,mean.July.12,mean.July.13,mean.July.14,mean.July.15,mean.July.16)

April.May.sst<-c(mean.April.May.78,mean.April.May.79,mean.April.May.80,mean.April.May.81,mean.April.May.82,mean.April.May.83,mean.April.May.84,mean.April.May.85,
                 mean.April.May.86,mean.April.May.87,mean.April.May.88,mean.April.May.89,mean.April.May.90,mean.April.May.91,mean.April.May.92,mean.April.May.93,mean.April.May.94,
                 mean.April.May.95,mean.April.May.96,mean.April.May.97,mean.April.May.98,mean.April.May.99,mean.April.May.00,mean.April.May.01,mean.April.May.02,mean.April.May.03,
                 mean.April.May.04,mean.April.May.05,mean.April.May.06,mean.April.May.07,mean.April.May.08,mean.April.May.09,mean.April.May.10,mean.April.May.11,mean.April.May.12,mean.April.May.13,
                 mean.April.May.14,mean.April.May.15,mean.April.May.16)

May.June.sst<-c(mean.May.June.78,mean.May.June.79,mean.May.June.80,mean.May.June.81,mean.May.June.82,mean.May.June.83,mean.May.June.84,mean.May.June.85,mean.May.June.86,
                mean.May.June.87,mean.May.June.88,mean.May.June.89,mean.May.June.90,mean.May.June.91,mean.May.June.92,mean.May.June.93,mean.May.June.94,mean.May.June.95,mean.May.June.96,
                mean.May.June.97,mean.May.June.98,mean.May.June.99,mean.May.June.00,mean.May.June.01,mean.May.June.02,mean.May.June.03,mean.May.June.04,mean.May.June.05,mean.May.June.06,
                mean.May.June.07,mean.May.June.08,mean.May.June.09,mean.May.June.10,mean.May.June.11,mean.May.June.12,mean.May.June.13,mean.May.June.14,mean.May.June.15,mean.May.June.16)


June.July.sst<-c(mean.June.July.78,mean.June.July.79,mean.June.July.80,mean.June.July.81,mean.June.July.82,mean.June.July.83,mean.June.July.84,mean.June.July.85,mean.June.July.86,
                 mean.June.July.87,mean.June.July.88,mean.June.July.89,mean.June.July.90,mean.June.July.91,mean.June.July.92,mean.June.July.93,mean.June.July.94,mean.June.July.95,mean.June.July.96,
                 mean.June.July.97,mean.June.July.98,mean.June.July.99,mean.June.July.00,mean.June.July.01,mean.June.July.02,mean.June.July.03,mean.June.July.04,mean.June.July.05,mean.June.July.06,
                 mean.June.July.07,mean.June.July.08,mean.June.July.09,mean.June.July.10,mean.June.July.11,mean.June.July.12,mean.June.July.13,mean.June.July.14,mean.June.July.15,mean.June.July.16)


April.July.sst<-c(mean.April.July.78,mean.April.July.79,mean.April.July.80,mean.April.July.81,mean.April.July.82,mean.April.July.83,mean.April.July.84,mean.April.July.85,
                  mean.April.July.86,mean.April.July.87,mean.April.July.88,mean.April.July.89,mean.April.July.90,mean.April.July.91,mean.April.July.92,mean.April.July.93,mean.April.July.94,
                  mean.April.July.95,mean.April.July.96,mean.April.July.97,mean.April.July.98,mean.April.July.99,mean.April.July.00,mean.April.July.01,mean.April.July.02,mean.April.July.03,
                  mean.April.July.04,mean.April.July.05,mean.April.July.06,mean.April.July.07,mean.April.July.08,mean.April.July.09,mean.April.July.10,mean.April.July.11,mean.April.July.12,mean.April.July.13,
                  mean.April.July.14,mean.April.July.15,mean.April.July.16)


May.July.sst<-c(mean.May.July.78,mean.May.July.79,mean.May.July.80,mean.May.July.81,mean.May.July.82,mean.May.July.83,mean.May.July.84,mean.May.July.85,mean.May.July.86,
                mean.May.July.87,mean.May.July.88,mean.May.July.89,mean.May.July.90,mean.May.July.91,mean.May.July.92,mean.May.July.93,mean.May.July.94,mean.May.July.95,mean.May.July.96,
                mean.May.July.97,mean.May.July.98,mean.May.July.99,mean.May.July.00,mean.May.July.01,mean.May.July.02,mean.May.July.03,mean.May.July.04,mean.May.July.05,mean.May.July.06,
                mean.May.July.07,mean.May.July.08,mean.May.July.09,mean.May.July.10,mean.May.July.11,mean.May.July.12,mean.May.July.13,mean.May.July.14,mean.May.July.15,mean.May.July.16)





#############################################################################################################################################
#################################################### Plot ###################################################################################
par(mfrow=c(2,2))
Year<-c(1978:2016)
plot(April.sst~Year,pch=16)
plot(May.sst~Year,pch=16)
plot(June.sst~Year,pch=16)
plot(July.sst~Year,pch=16)
plot(April.May.sst~Year,pch=16)
plot(May.June.sst~Year,pch=16)
plot(June.July.sst~Year,pch=16)
plot(April.July.sst~Year,pch=16)
plot(May.July.sst~Year,pch=16)



#############################################################################################################################################
#################################################### Calculate anomolies ####################################################################
#############################################################################################################################################
April.series.mean<-mean(April.sst)
April.series.mean
April.anom<-April.sst-April.series.mean
April.sst
April.anom

May.series.mean<-mean(May.sst)
May.series.mean
May.anom<-May.sst-May.series.mean
May.sst
May.anom

June.series.mean<-mean(June.sst)
June.series.mean
June.anom<-June.sst-June.series.mean
June.sst
June.anom

July.series.mean<-mean(July.sst)
July.series.mean
July.anom<-July.sst-July.series.mean
July.sst
July.anom

April.May.series.mean<-mean(April.May.sst)
April.May.series.mean
April.May.anom<-April.May.sst-April.May.series.mean
April.May.sst
April.May.anom

May.June.series.mean<-mean(May.June.sst)
May.June.series.mean
May.June.anom<-May.June.sst-May.June.series.mean
May.June.sst
May.June.anom

June.July.series.mean<-mean(June.July.sst)
June.July.series.mean
June.July.anom<-June.July.sst-June.July.series.mean
June.July.sst
June.July.anom

April.July.series.mean<-mean(April.July.sst)
April.July.series.mean
April.July.anom<-April.July.sst-April.July.series.mean
April.July.sst
April.July.anom

May.July.series.mean<-mean(May.July.sst)
May.July.series.mean
May.July.anom<-May.July.sst-May.July.series.mean
May.July.sst
May.July.anom

#####################################################################################################################################
############################################### Plot anomalies ######################################################################
plot(April.anom, pch=16)
plot(May.anom, pch=16)
plot(June.anom, pch=16)
plot(July.anom, pch=16)
plot(April.May.anom, pch=16)
plot(May.June.anom, pch=16)
plot(June.July.anom, pch=16)
plot(April.July.anom, pch=16)
plot(May.July.anom, pch=16)

####################################### Series for multivasriate analyses ######################################
################################################################################################################
April.anom
May.anom
June.anom
July.anom
April.May.anom
May.June.anom
June.July.anom
April.July.anom
May.July.anom
#################################################################################################################

sst.dat<-cbind(Year, April.sst, May.sst, June.sst, July.sst, April.May.sst, May.June.sst, June.July.sst, April.July.sst, May.July.sst)

sst.anom.dat<-cbind(Year,April.anom, May.anom, June.anom, July.anom, April.May.anom, May.June.anom, June.July.anom, April.July.anom, May.July.anom)

sst.comb.dat<-cbind(Year, April.sst, May.sst, June.sst, July.sst, April.May.sst, May.June.sst, June.July.sst, April.July.sst, May.July.sst,
April.anom, May.anom, June.anom, July.anom, April.May.anom, May.June.anom, June.July.anom, April.July.anom, May.July.anom)


write.csv(sst.comb.dat,"data/ERSST_SST_avgs_anoms.csv")

                  

