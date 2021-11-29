#setwd("C:\\Users\\User\\Dropbox\\CLOexternalBrain\\projects\\MA Thesis\\code\\data")
setwd("/Users/christoffer/Dropbox/CLOexternalBrain/projects/MA Thesis/code/data")
library(data.table)
library(ggplot2)
library(ggpubr)

states_hard = fread("raw/state_data_hard.csv")
states_easy = fread("raw/state_data_easy.csv")

stream_hard = fread("results/causal_streams_hard.csv")
stream_easy = fread("results/causal_streams_easy.csv")

stream_plot = function(states, stream, trial, start, stop, limit, title_label, no_stream = F, trials = c(6,24,38,57,70,91,102,122)){
  
  t = trials[trial]
  time = limit:stop
  states_sub = states[trial == t & timestep %in% time, 2:10]
  
  t = trial
  stream_sub = stream[trial == t & persp_element %in% 3:4 & persp_timestep %in% start:stop & cause_timestep >= limit]
  
  from = c(1,2,3,4,5,6,7,8)
  to = c(8,7,2,1,6,5,4,3)
  
  stream_sub[ , reader := to[match(reader, from)]]
  stream_sub[ , cause := to[match(cause, from)]]
  
  stream_sub[ , persp_element := c("M1","M2")[match(persp_element, c(3,4))]]
  
  half_range = (stop-start)/10/2
  stream_sub$cause = stream_sub$cause - half_range
  stream_sub$reader = stream_sub$reader - half_range
  stream_sub$cause_timestep = stream_sub$cause_timestep - half_range
  stream_sub$reader_timestep = stream_sub$reader_timestep - half_range
  
  point_size = 10
  
  element_labels = data.frame(
    x = limit-0.75,
    y = to,
    labels = c("S1", "S2", "M1", "M2", "H1", "H2", "H3", "H4")
  )
  
  plot = ggplot(stream_sub) +
    geom_point(data = states_sub, aes(x=timestep, y = 8, fill = as.factor(S1)), shape = 21, size=point_size) +
    geom_point(data = states_sub, aes(x=timestep, y = 7, fill = as.factor(S2)), shape = 21, size=point_size) +
    geom_point(data = states_sub, aes(x=timestep, y = 2, fill = as.factor(M1)), shape = 21, size=point_size) +
    geom_point(data = states_sub, aes(x=timestep, y = 1, fill = as.factor(M2)), shape = 21, size=point_size) +
    geom_point(data = states_sub, aes(x=timestep, y = 6, fill = as.factor(H1)), shape = 21, size=point_size) +
    geom_point(data = states_sub, aes(x=timestep, y = 5, fill = as.factor(H2)), shape = 21, size=point_size) +
    geom_point(data = states_sub, aes(x=timestep, y = 4, fill = as.factor(H3)), shape = 21, size=point_size) +
    geom_point(data = states_sub, aes(x=timestep, y = 3, fill = as.factor(H4)), shape = 21, size=point_size) +
    scale_fill_manual(values = c("white","black")) +
    
    theme_minimal() +
    theme(legend.position = "none",
          panel.spacing = unit(3, "lines"),
          strip.text.x = element_text(size = 20, face = "bold"),
          panel.grid.minor.x = element_line(color = "gray",
                                            size = 0.5,
                                            linetype = 2),
          panel.grid.major.x = element_line(color = "white")
          ) +
    labs(title = title_label, y = " ", x = " ")+
    scale_x_continuous(breaks = time) +
    scale_y_continuous(breaks=NULL) +
    geom_text(data = element_labels, aes(x = x, y = y, label = labels), size=4)+
    geom_segment(data = element_labels, aes(x=limit-1, xend=stop+.5, y=y-.5, yend=y-.5), color = "lightgray")
  
  if(no_stream == F){
    plot = plot + 
      facet_wrap(~persp_element) +
      geom_segment(aes(y = cause+(persp_timestep-start)/10, 
                       yend = reader+(persp_timestep-start)/10, 
                       x = cause_timestep+(persp_timestep-start)/10, 
                       xend = reader_timestep+(persp_timestep-start)/10, 
                       color = as.factor(persp_timestep)),
                   size = 1, alpha = .7) +
      geom_point(aes(y = cause+(persp_timestep-start)/10, 
                     x = cause_timestep+(persp_timestep-start)/10, 
                     color = as.factor(persp_timestep)),
                 size = 5)+
      geom_point(aes(y = reader+(persp_timestep-start)/10, 
                     x = reader_timestep+(persp_timestep-start)/10, 
                     color = as.factor(persp_timestep)),
                 size = 5)
  }
  
  
  return(plot)
  
}

stream_plot(states_easy, stream_easy, 4, 12, 15, 3, "EASY TASK\nBlock size: 3\nTrial type: Avoid\nDirection: Left", F, c(6,24,39,57))

ggsave("../figures/trial_states_easy.jpg",
  ggarrange(
    stream_plot(states_easy, stream_easy, 3, 1, 15, 5, "EASY TASK\nBlock size: 3\nTrial type: Avoid\nDirection: Left", T, c(6,24,39,57)),
    stream_plot(states_easy, stream_easy, 4, 1, 15, 5, "EASY TASK\nBlock size: 3\nTrial type: Avoid\nDirection: Right", T, c(6,24,39,57)),
    nrow = 2, ncol = 1
  ), height = 10, width = 8
)

ggsave("../figures/trial_states_hard.jpg",
       ggarrange(
         stream_plot(states_hard, stream_hard, 1, 1, 18, 5, "HARD TASK\nBlock size: 3\nTrial type: Catch\nDirection: Left", T),
         stream_plot(states_hard, stream_hard, 4, 1, 18, 5, "HARD TASK\nBlock size: 4\nTrial type: Avoid\nDirection: Right", T),
         nrow = 2, ncol = 1
       ), height = 10, width = 10
)

ggsave("../figures/streams_easy1.jpg",
       stream_plot(states_easy, stream_easy, 3, 8, 12, 3, "EASY TASK\nBlock size: 3\nTrial type: Avoid\nDirection: Left", F, c(6,24,39,57)),
       height = 6, width = 15
)

ggsave("../figures/streams_easy2.jpg",
       stream_plot(states_easy, stream_easy, 4, 9, 13, 3, "EASY TASK\nBlock size: 3\nTrial type: Avoid\nDirection: Right", F, c(6,24,39,57)),
       height = 6, width = 15
)

ggsave("../figures/streams_hard1.jpg",
       stream_plot(states_hard, stream_hard, 1, 10, 13, 6, "HARD TASK\nBlock size: 3\nTrial type: Catch\nDirection: Left"),
       height = 6, width = 15
)

ggsave("../figures/streams_hard2.jpg",
       stream_plot(states_hard, stream_hard, 4, 9, 12, 5, "HARD TASK\nBlock size: 4\nTrial type: Avoid\nDirection: Right"),
       height = 6, width = 15
)
