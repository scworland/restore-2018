
set.seed(1)
centers <- data.frame(cluster=factor(1:3), 
                      size=c(100, 150, 50), 
                      x1=c(5, 0, -3), 
                      x2=c(-1, 1, -2))

year1 <- centers %>% 
  group_by(cluster) %>%
  do(data.frame(x1=rnorm(.$size[1], .$x1[1]),
                x2=rnorm(.$size[1], .$x2[1]),
                year="year 1",
                stringsAsFactors = F)) %>%
  data.frame()

year2 <- centers %>% 
  group_by(cluster) %>%
  do(data.frame(x1=rnorm(.$size[1], .$x1[1]),
                x2=rnorm(.$size[1], .$x2[1]),
                year="year 2",
                stringsAsFactors = F)) %>%
  data.frame()


points <- rbind(year1,year2)

kclusters <- points %>%
  select(-cluster) %>%
  group_by(year) %>%
  do(data.frame(., kclust = kmeans(as.matrix(.[,-3]),centers=3)$cluster)) %>%
  mutate(kclust = as.character(kclust))

ggplot(kclusters) + 
  geom_point(aes(x1,x2,color=kclust)) + 
  facet_wrap(~year) +
  theme_bw() +
  scale_color_viridis_d()
