library(EBImage)
library(matlab)
library(ncdf4)
library(maps)
library(leaflet)
rm(list = ls())




lx = seq(90,150,length=3002)
ly = seq(-15,15,length=1502)
ilx = which(lx >= 115 & lx <= 125)
ily = which(ly >= -7 & ly <= 0)
setwd("~/Desktop/Project_BW4/program_10mnt_mtsat/OUTPUT/")
load("~/Desktop/Studi_Yosik/data/topograp.rda")
# setwd("~/Desktop/Studi_Yosik/r_code/proven_method/grab-tag-graph-master/datadir/MERG")
f = list.files(pattern = "nc")
data_nc=f[length(f)]
# batas_luas = 1000
pilih_ecc  =0.7

filter_A_ECC = function(data_nc,batas_luas,pilih_ecc_min,pilih_ecc_max,ilx,ily){
  nc = nc_open(data_nc)
  var = ncvar_get(nc,"suhu")[ilx,ily]
  var[var >= 250]=0
  lon = ncvar_get(nc,"lon")[ilx]
  lat = ncvar_get(nc,"lat")[ily]
  length(lat)*length(lon)
  dis = lat[length(lat)] - lat[1]
  degree = dis/length(lat)
  const = 111.699
  km2 = batas_luas
  jumlah_grid = as.integer(km2/(degree*const))
  
  l<-bwlabel(var)
  idI=l
  l[l == 0] = NA
  
  
  shape_var = computeFeatures.shape(var)
  shape_var = as.data.frame(shape_var)
  n_shape = names(shape_var)
  moment_var = computeFeatures.moment(l)
  moment_var = as.data.frame(moment_var)
  n_moment = names(moment_var)
  row_shape = rownames(shape_var)
  row_mom = rownames(moment_var)
  
  RegionProps <- table(l)
  length(RegionProps)
  names(RegionProps)
  
  # pilih_ecc = 0.7
  iecc = which(moment_var$m.eccentricity > pilih_ecc_min & moment_var$m.eccentricity <= pilih_ecc_max)   #yang bulat dibuang
  simpan_ecc =  !which(moment_var$m.eccentricity > pilih_ecc_min & moment_var$m.eccentricity <= pilih_ecc_max )
  AreaThreshold = jumlah_grid # jumlah grid
  idx<-which(RegionProps < AreaThreshold) #delete those
  simpan_luas<-row_shape[shape_var$s.area >= AreaThreshold] #delete those
  
  simpan_ecc1 = simpan_ecc[simpan_ecc %in% as.numeric(simpan_luas)]
  simpan_luas1 = simpan_luas[simpan_luas %in% simpan_ecc1]
  simpan_luas1 =  shape_var[row_shape %in% simpan_luas1,]
  simpan_ecc1 = moment_var[row_mom %in% simpan_ecc1,]
  hasil_doc = data.frame(simpan_luas1,simpan_ecc1)
  
  
  retImage<-var
  
  for (i in iecc){
    # print(paste(i,' of ', length(iecc),sep=''))
    retImage[l == i ]<- 0
    idI[l == i ]<- 0
  }
  for (i in idx){
    # print(paste(i,' of ', length(iecc),sep=''))
    retImage[l == i ]<- 0
    idI[l == i ] <- 0
  }
  data_doc = list(retImage,idI)
  
  doc_area <- RegionProps[RegionProps > AreaThreshold]
  
  doc_area = doc_area[names(doc_area) != "0"]
  
  
  save_doc = moment_var[which(rownames(moment_var) %in% names(doc_area)),]
  doc_area = as.numeric(doc_area)
  save_doc = cbind(doc_area,save_doc)
  
  
  unique(retImage)
  names(data_doc) = c("brightnesstemp","cluster")
  hasil = list(data_doc,save_doc)
  # hasil = list(data_doc,hasil_doc)
  names(hasil) = c("data","documentation")
  xx = hasil$documentation$m.cx
  yy = hasil$documentation$m.cy
  # hasil$documentation$m.cx = ((xx/length(lon))*(lon[length(lon)] - lon[1]))-lon[length(lon)]
  # hasil$documentation$m.cy = ((xx/length(lat))*(lat[length(lat)] - lat[1]))+lat[1]
  return(hasil)
}

function(data_nc,batas_luas,pilih_ecc_min,pilih_ecc_max,ilx,ily,core_center_threshold  ){
  oke = filter_A_ECC(data_nc = data_nc,batas_luas = batas_luas,pilih_ecc_min = pilih_ecc_min ,
                     pilih_ecc_max = pilih_ecc_max,ilx,ily)
  num_dat = as.numeric(oke$data[[2]])  
  i_x=unique(num_dat)
  bt = oke$data[[1]]
  xoxo = bt
  i_bt = oke$data[[2]]
  cory = ly[ily]
  corx = lx[ilx]
  xy = meshgrid(1:length(cory),1:length(corx))
  
  xoxo = list()
  min_im = c()
  c_x = c()
  c_y = c()
  yang_kepilih = c()
  for(i in 1:length(i_x)){
    xoxo[[i]]=bt
    xoxo[[i]][i_bt != i_x[i]] = 0
    min_im[i] = min(xoxo[[i]][xoxo[[i]] !=0])
    if(min_im[i] <= 250){
      c_x[i] = xy$y[xoxo[[i]] == min_im[i]]
      c_y[i] = xy$x[xoxo[[i]] == min_im[i]]
    }
  }
  rangeee= cory[length(cory)]-cory[1]
  terpilih = which(!is.na(c_x))
  # uniq_code = 
  
  hasil = list()
  hasil$doc = oke$documentation[terpilih,]
  
  c_y1 = cory[c_y]
  c_x1 = corx[c_x]
  xxyy = bt
  xxyy[bt ==0]=NA
  hasil$doc$m.cx = c_x1[!is.na(c_x1)]
  hasil$doc$m.cy = c_y1[!is.na(c_y1)]
  hasil$bt = xxyy
  hasil$lon = corx
  hasil$lat = cory
  return(hasil)
}


oke = filter_A_ECC(data_nc = f[9],batas_luas = 10000,pilih_ecc_min = 0.7,pilih_ecc_max = 1,ilx,ily)

num_dat = as.numeric(oke$data[[2]])  
i_x=unique(num_dat)
bt = oke$data[[1]]
xoxo = bt

i_bt = oke$data[[2]]
cory = ly[ily]
corx = lx[ilx]
xy = meshgrid(1:length(cory),1:length(corx))
length(corx)


xoxo = list()
min_im = c()
c_x = c()
c_y = c()
yang_kepilih = c()
for(i in 1:length(i_x)){
  xoxo[[i]]=bt
  xoxo[[i]][i_bt != i_x[i]] = 0
  min_im[i] = min(xoxo[[i]][xoxo[[i]] !=0])
  if(min_im[i] <= 250){
    c_x[i] = xy$y[xoxo[[i]] == min_im[i]]
    c_y[i] = xy$x[xoxo[[i]] == min_im[i]]
  }
}
# image(1:414,1:414,bt)


rangeee= cory[length(cory)]-cory[1]

# points(c_x,c_y)
terpilih = which(!is.na(c_x))
# uniq_code = 

hasil = list()
hasil$doc = oke$documentation[terpilih,]

c_y1 = cory[c_y]
c_x1 = corx[c_x]
xxyy = bt
xxyy[bt ==0]=NA

# image(xnya,ynya,el,axes=FALSE,xlab='Longitude',ylab='Latitude',col=grey.colors(100))

colors = c(rgb(0.5,0.2,0.7,alpha = 0.7),
           rgb(1,0.2,0.2,alpha = 0.7), 
           rgb(1,0.7,0.2,alpha = 0.7), 
           rgb(0.4,0.7,0.2,alpha = 0.7),
           rgb(1,0.2,0.9,alpha = 0.7),
           rgb(0.1,0.2,0.9,alpha = 0.7))
ccol = colorRampPalette(colors = c("white","grey","darkblue","orange","darkred","grey"))

x11()

filled.contour(corx,cory,xxyy,col = rev(ccol(50)),levels = seq(150,250,length=40),axes = T,
               plot.axes = {axis(1,cex.axis = 1.2);axis(2,cex.axis = 1.2);
                 title(main=sprintf("Convective Perimeter %s",substr(f[length(f)],1,12)),cex.main=3);
                 image(topograp$xnya,topograp$ynya,topograp$el,col=terrain.colors(10,alpha = 0.01),add = T);
                 grid(nx = 20,ny=15);points(c_x1,c_y1,col="black",pch=7,cex=2,lwd=2);
                 map("world", fill=F, col="brown", bg=NULL, xlim=c(90,150), ylim=c(-15, 15), 
                     mar=c(0,0,0,0),resolution = 0.0000001,add=T)})
CL <- contourLines(corx , cory , xxyy)
library(sp)
pgons <- lapply(1:length(CL), function(i)
  Polygons(list(Polygon(cbind(CL[[i]]$x, CL[[i]]$y))), ID=i))
spgons = SpatialPolygons(pgons)
LEVS <- as.factor(sapply(CL, `[[`, "level"))
NLEV <- length(levels(LEVS))
m=leaflet()
m = addTiles(m)
m = addMarkers(m, lng=c_x1[2:3], lat=c_y1[2:3], popup="ajozz")
leaflet(spgons) %>% addTiles() %>% 
  addPolygons(color = heat.colors(NLEV, NULL)[LEVS])  %>%
  addMarkers(m, lng=c_x1[2:3], lat=c_y1[2:3], popup=sprintf("Monitoring Kluster Awan Konvektif; data tanggal %s; jam %s UTC; pusat lokasi lon,lat : (%s,%s) ;
                                                            Luasan : %s Km^2; Perimeter : %s",
                                                            substr(data_nc,1,8),substr(data_nc,9,12),
                                                            substr(rev(c_x1[!is.na(c_x1)]),1,6),
                                                            substr(rev(c_y1[!is.na(c_y1)]),1,6),
                                                            rev(hasil$doc$doc_area),
                                                            substr(rev(hasil$doc$m.eccentricity),1,6)
                                                          )
             )
