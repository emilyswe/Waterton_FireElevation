packs<-c("runjags","beepr","BiocManager")
for(q in 1:length(packs)) 
	if (!require(packs[q], character.only = TRUE)) 
		{
			install.packages(packs[q])  
			require(packs[q])
		}
		
if (!require("ComplexHeatmap", character.only = TRUE)) BiocManager::install("ComplexHeatmap")
		
csv<-read.csv("04_Cleaned_Waterton.csv")

colu<-colnames(csv)
fsp<-which(colu=="COLO")
lsp<-which(colu=="STGR")
sps<-sort(colu[fsp:lsp])
wsps<-1:length(sps)
hml<-length(unique(csv$location))
hms<-1+(lsp-fsp)
hmo<-length(csv[,lsp])
hmy<-max(csv$year-2006)

###################### one model, species as random effect
####### get land cover array
dic<-sy<-lanc<-lc<-ele<-slo<-sev<-wsp<-asp<-tpi<-changeflag<-resparray<-north<-yr<-rep(0,2)


locs<-sort(unique(csv$location))

for(q in 1:hml) sy[q]<-min(csv[which(csv$location==locs[q]),which(colu=="year")])

for(q in 1:hmo){
	print(which(colu==paste0("landCover_",csv$year[q])))
	wc<-ifelse(length(which(colu==paste0("landCover_",csv$year[q])))>0,which(colu==paste0("landCover_",csv$year[q])),145)
	lc[q]<-csv[q,wc]
}

hlc<-matrix(ncol=hmy,nrow=hml,dimnames=list(locs,2007:2022),data=NA)
for(q in 1:hml)
	for(w in 1:(sy[q]-2006)) hlc[q,w]<-0
hlc[,1]<-rep(0,hml)
for(q in 1:hml)
	for(w in (sy[q]-2005):hmy)
		hlc[q,w]<-ifelse(length(csv[which(csv$location==rownames(hlc)[q]&csv$year==(w+2006)),145])>0,csv[which(csv$location==rownames(hlc)[q]&csv$year==(w+2006)),145],NA)

for(q in 1:hml)
	for(w in 2:13)
		if(is.na(hlc[q,w])==TRUE) hlc[q,w]<-hlc[q,w-1]

hlw<-hlc

for(q in 1:hml)
	for(w in 2:13)
		hlw[q,w]<-ifelse(hlc[q,w]==hlc[q,w-1],0,1)
		

for(q in 1:hmo)
	dic[q]<-hlw[which(rownames(hlc)==csv$location[q]),which(colnames(hlw)==csv$year[q])]

i<-1
for(q in 1:hms)
	{
	resparray[i:(i+hmo-1)]<-csv[,(q+fsp-1)]
	print(length(resparray))
	wsp[i:(i+hmo-1)]<-rep(which(sps==colu[(q+fsp-1)]),hmo)
	print(length(wsp))
	yr[i:(i+hmo-1)]<-csv$year-2006
	print(length(yr))
	lanc[i:(i+hmo-1)]<-lc
	print(length(lanc))
	sev[i:(i+hmo-1)]<-csv$gridcode+1
	print(length(sev))
	ele[i:(i+hmo-1)]<-scale(csv$elevation)[,1]
	print(length(ele))
	slo[i:(i+hmo-1)]<-scale(csv$slope)[,1]
	print(length(slo))
	asp[i:(i+hmo-1)]<-scale(csv$aspect)[,1]
	print(length(asp))
	tpi[i:(i+hmo-1)]<-scale(csv$TPI)[,1]
	print(length(tpi))
	changeflag[i:(i+hmo-1)]<-dic
	print(length(changeflag))
	north[i:(i+hmo-1)]<-scale(csv$northness)[,1]
	print(length(north))
	print(q)
#	print(colu[(q+fsp-1)])
	i<-i+hmo
	}
	
nas<-which(is.na(lanc))
resparray<-resparray[-nas]
wsp<-wsp[-nas]
yr<-yr[-nas]
sev<-sev[-nas]
ele<-ele[-nas]
slo<-slo[-nas]
asp<-asp[-nas]
tpi<-tpi[-nas]
lanc<-lanc[-nas]
changeflag<-changeflag[-nas]
north<-north[-nas]

	print(length(resparray))
		print(length(wsp))
			print(length(yr))
				print(length(lanc))
					print(length(sev))
						print(length(ele))
							print(length(slo))
								print(length(asp))
									print(length(tpi))
										print(length(changeflag))
											print(length(north))
											
lanc_types<-sort(unique(lanc))
nhmo<-length(asp)
for(q in 1:nhmo) lanc[q]<-which(lanc_types==lanc[q])
datalist<-list(resp=resparray,wsp=wsp,nsp=hms,hmo=nhmo,yr=yr,ny=max(yr),nlc=length(lanc_types),lc=lanc,nsev=max(sev),sev=sev,ele=ele,slo=slo,asp=asp,tpi=tpi,cf=changeflag+1,north=north)


# Dummy variables for JAGS settings
n_chains <- 5 # Number of chains
n_iter <- 5e4 # Number of iterations
n_thin <- 5 # Thinning rate

n_cores <- n_chains 

keep_track_of <- c("alpha","alpha_sp","beta_yr","beta_lc","beta_sev","beta_cf","beta_ele","beta_slo","beta_asp","beta_tpi","beta_north")
  
model_file='rmod.txt'
# Run the model with runjags
result <- runjags::run.jags(
	  model = model_file,
	  monitor = keep_track_of,
	  data = datalist,
	  n.chains = n_chains,
	  adapt = n_iter/10,
	  burnin = n_iter*2, 
	  sample = n_iter/(n_thin), 
	  thin = n_thin,
	  method = "parallel"
	)
save.image(file="rmod.Rdata",version=2)
# Print the summary of the result
print(result)
sure<-summary(result)

suva<-rownames(sure)

lanc_names<-rep(0,2)
for(q in 1:length(lanc_types)) lanc_names[q]<-csv$category[which(csv$code==lanc_types[q])][1]
heatcolumns<-c("alpha",2007:2019,lanc_names,paste("burn severity",0:5),"elevation","slope","aspect","TPI","Northness","Same Land Cover","Change in Land Cover")
hmhc<-length(heatcolumns)
reta<-matrix(nrow=length(wsps),ncol=hmhc,dimnames=list(sps,heatcolumns))

alphas<-grep("alpha_sp",suva)
reta[,1]<-sure[alphas,2]

cv<-grep("beta_yr",suva)

for(q in 1:hms)
	for(w in 1:max(yr))
		reta[q,1+w]<-sure[(cv[q]+111*(w-1)),2]

cv<-grep("beta_lc",suva)

for(q in 1:hms)
	for(w in 1:length(lanc_names))
		reta[q,14+w]<-sure[(cv[q]+111*(w-1)),2]

cv<-grep("beta_sev",suva)

for(q in 1:hms)
	for(w in 1:length(lanc_names))
		reta[q,21+w]<-sure[(cv[q]+111*(w-1)),2]

#28
w<-28
cv<-grep("beta_ele",suva)
reta[,w]<-sure[cv,2]

w<-w+1
cv<-grep("beta_slo",suva)
reta[,w]<-sure[cv,2]

w<-w+1
cv<-grep("beta_asp",suva)
reta[,w]<-sure[cv,2]

w<-w+1
cv<-grep("beta_tpi",suva)
reta[,w]<-sure[cv,2]

w<-w+1
cv<-grep("beta_north",suva)
reta[,w]<-sure[cv,2]


cv<-grep("beta_cf",suva)
for(q in 1:hms)
	for(e in 1:2)
		reta[q,e+w]<-sure[(cv[q]+111*(e-1)),2]


###########


# Initialize the binary_matrix with dimensions and names identical to 'reta'
binary_matrix <- matrix(nrow=length(wsps), ncol=hmhc, dimnames=list(sps, heatcolumns))

# Process 'alpha' values
alphas <- grep("alpha_sp", suva)
binary_matrix[,1] <- ifelse(sure[alphas,1] > 0 | sure[alphas,3] < 0, 1, 0)

# Process 'beta_yr' values
cv <- grep("beta_yr", suva)
for(q in 1:hms) {
  for(w in 1:max(yr)) {
    idx <- cv[q] + 111 * (w - 1)
    binary_matrix[q, 1 + w] <- ifelse(sure[idx, 1] > 0 | sure[idx, 3] < 0, 1, 0)
  }
}

# Process 'beta_lc' values for land cover categories
cv <- grep("beta_lc", suva)
for(q in 1:hms) {
  for(w in 1:length(lanc_names)) {
    idx <- cv[q] + 111 * (w - 1)
    binary_matrix[q, 14 + w] <- ifelse(sure[idx, 1] > 0 | sure[idx, 3] < 0, 1, 0)
  }
}

# Process 'beta_sev' values for burn severity categories
cv <- grep("beta_sev", suva)
for(q in 1:hms) {
  for(w in 1:6) {  # Assuming there are 6 burn severity categories
    idx <- cv[q] + 111 * (w - 1)
    binary_matrix[q, 21 + w] <- ifelse(sure[idx, 1] > 0 | sure[idx, 3] < 0, 1, 0)
  }
}

# Variables 'elevation', 'slope', 'aspect', 'TPI', 'Northness'
# For each of these variables, there's only one corresponding beta value per species

# Process 'beta_ele' (elevation)
w <- 28  # Assuming 'elevation' is column 28 based on your setup
cv <- grep("beta_ele", suva)
binary_matrix[, w] <- ifelse(sure[cv, 1] > 0 | sure[cv, 3] < 0, 1, 0)

# Process 'beta_slo' (slope)
w <- w + 1  # Update column index for 'slope'
cv <- grep("beta_slo", suva)
binary_matrix[, w] <- ifelse(sure[cv, 1] > 0 | sure[cv, 3] < 0, 1, 0)

# Process 'beta_asp' (aspect)
w <- w + 1  # Update column index for 'aspect'
cv <- grep("beta_asp", suva)
binary_matrix[, w] <- ifelse(sure[cv, 1] > 0 | sure[cv, 3] < 0, 1, 0)

# Process 'beta_tpi' (TPI)
w <- w + 1  # Update column index for 'TPI'
cv <- grep("beta_tpi", suva)
binary_matrix[, w] <- ifelse(sure[cv, 1] > 0 | sure[cv, 3] < 0, 1, 0)

# Process 'beta_north' (Northness)
w <- w + 1  # Update column index for 'Northness'
cv <- grep("beta_north", suva)
binary_matrix[, w] <- ifelse(sure[cv, 1] > 0 | sure[cv, 3] < 0, 1, 0)

# Process 'Same Land Cover' and 'Change in Land Cover' conditions
# Assuming these correspond to the next two columns after 'Northness'

# Process 'beta_cf' for 'Same Land Cover'
w <- w + 1  # Update column index for 'Same Land Cover'
cv <- grep("beta_cf", suva)
for(q in 1:hms) {
    idx <- cv[q]  # Assuming only one value for 'Same Land Cover'
    binary_matrix[q, w] <- ifelse(sure[idx, 1] > 0 | sure[idx, 3] < 0, 1, 0)
}

# Increment w for 'Change in Land Cover'
w <- w + 1
for(q in 1:hms) {
    idx <- cv[q] + 111  # Assuming the next index for 'Change in Land Cover'
    binary_matrix[q, w] <- ifelse(sure[idx, 1] > 0 | sure[idx, 3] < 0, 1, 0)
}

lbrary(ComplexHeatmap)

jm<- Heatmap(reeta, 
        cluster_columns = FALSE, 
        cluster_rows = FALSE, 
        show_column_names = TRUE, 
        show_row_names = TRUE, 
        heatmap_legend_param = list(title = "Median Value", title_position = "topcenter"))











##################### one model: multinomial


