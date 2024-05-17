theme_adlunglab <- function(base_size = 12, base_family = "Helvetica", 
                            legend_position = "right", title_size= base_size){
theme_bw(base_size = base_size, base_family = base_family) %+replace%
theme(legend.position = legend_position, legend.background = element_blank(),
strip.background = element_blank(),
panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
axis.line = element_line(colour = "#000000"), axis.ticks = element_line(colour = "#000000"),
legend.key = element_blank(),
axis.text = element_text(size = base_size), plot.title=element_text(size = title_size),
axis.title = element_text(size = base_size), legend.text=element_text(size = base_size),
legend.title=element_text(size = base_size), strip.text=element_text(size = base_size)
)
}


ALcols <- c("#D8D495", "#D4D958", "#7BC086", "#09A64A", "#00B2D9", "#8D2668", "#8C1912", "#D90767", "#EB5D12", "#F39F07")
