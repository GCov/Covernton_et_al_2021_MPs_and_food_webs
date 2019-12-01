sea_cucumbers <- read.csv('sea_cucumbers.csv')

sea_cucumbers$log <- with(sea_cucumbers,
                          ifelse(colour == 'blue' &
                                   shape == 'fragment',
                                 'TRUE', 'FALSE'))

sea_cucumbers2 <- subset(sea_cucumbers, log != 'TRUE')

ID <- levels(sea_cucumbers$ID)

CUIDs <- data.frame(ID)

sea_cucumbers2 <- left_join(CUIDs, sea_cucumbers2, by = 'ID')

write.csv(sea_cucumbers2, 'sea_cucumbers.csv')
