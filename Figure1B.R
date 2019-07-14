Race_Distribution = sapply(split(qq_df$race, qq_df$hist), function(x) table(x)[table(x)>0])
df=melt(Race_Distribution, id=c("AA"))
colnames(df)=c('Race', 'Samples_Count', 'Cancer_Type')


tiff('/cbcb/project2-scratch/sanju/stack_delthis.tiff', height=800, width=1000)
ggplot(df, aes(Cancer_Type, y=Samples_Count))  +
	geom_bar(stat='identity', aes(fill=Race))+
	theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=25))
dev.off()


