library(fossil)
library(vegan)
library(dplyr)
library(reshape2)
#Geo. distance
locations <- read.table("D:/OneDrive - zju.edu.cn/Gym/environment/get_env_factor/gym_22pop_location.txt",header = TRUE) %>%tibble::column_to_rownames('pop') 
Geodistancematrix <- as.matrix( earth.dist(locations))
Fst_scalar<- read.table("D:/OneDrive - zju.edu.cn/Gym/Mantel test/gym_22pop_fstlongitude_scalar.txt",header = TRUE) %>% as.matrix
#Env. distance
all104env<- read.table("D:/OneDrive - zju.edu.cn/Gym/Mantel test/envfactor104_pop22_PCA_matrix.txt",header = TRUE)%>% as.matrix
longitude<- read.table("D:/OneDrive - zju.edu.cn/Gym/Mantel test/longitude_pop22_matrix.txt",header = TRUE)%>% as.matrix
latitude<- read.table("D:/OneDrive - zju.edu.cn/Gym/Mantel test/latitude_pop22_matrix.txt",header = TRUE)%>% as.matrix
pre<- read.table("D:/OneDrive - zju.edu.cn/Gym/Mantel test/all22pop_pre_PCA_matrix.txt",header = TRUE)%>% as.matrix
tem<- read.table("D:/OneDrive - zju.edu.cn/Gym/Mantel test/all22pop_tem_PCA_matrix.txt",header = TRUE)%>% as.matrix
srad<- read.table("D:/OneDrive - zju.edu.cn/Gym/Mantel test/all22pop_srad_PCA_matrix.txt",header = TRUE)%>% as.matrix
wind<- read.table("D:/OneDrive - zju.edu.cn/Gym/Mantel test/all22pop_wind_PCA_matrix.txt",header = TRUE)%>% as.matrix
vapr<- read.table("D:/OneDrive - zju.edu.cn/Gym/Mantel test/all22pop_vapr_PCA_matrix.txt",header = TRUE)%>% as.matrix
salt<- read.table("D:/OneDrive - zju.edu.cn/Gym/environment/get_env_factor/salt_matrix.txt",header = TRUE)%>% as.matrix
AI<- read.table("D:/OneDrive - zju.edu.cn/Gym/environment/get_env_factor/AI_matrix.txt",header = TRUE)%>% as.matrix
## Mantel test
mantel(Geodistancematrix, Fst_scalar)
mantel(Geodistancematrix, all104env) 
mantel(all104env, Fst_scalar)
mantel(longitude, Fst_scalar)
mantel(latitude, Fst_scalar)
mantel(pre, Fst_scalar)
mantel(tem, Fst_scalar)
mantel(srad, Fst_scalar)
mantel(vapr, Fst_scalar)
mantel(wind, Fst_scalar)
mantel(salt, Fst_scalar)
mantel(AI, Fst_scalar)
## partial Mantel test
mantel.partial(Geodistancematrix,Fst_scalar,all104env)
mantel.partial(all104env,Fst_scalar,Geodistancematrix)
mantel.partial(longitude, Fst_scalar,Geodistancematrix)
mantel.partial(latitude, Fst_scalar,Geodistancematrix)
mantel.partial(pre, Fst_scalar,Geodistancematrix)
mantel.partial(tem, Fst_scalar,Geodistancematrix)
mantel.partial(srad, Fst_scalar,Geodistancematrix)
mantel.partial(vapr, Fst_scalar,Geodistancematrix)
mantel.partial(wind, Fst_scalar,Geodistancematrix)
mantel.partial(salt, Fst_scalar,Geodistancematrix)
mantel.partial(AI, Fst_scalar,Geodistancematrix)