# barplot -----
bar_plot = function (object, probs, plot.flag,title, ...) {

  if (is.null(probs)) {
    object_summary = data.frame(mean = apply(object,2,mean, na.rm=T)  )
    object_summary$names = rownames(object_summary)
  } else if (is.null(probs) == F){
    object_summary = data.frame(mean = apply(object,2,mean, na.rm=T)   ,
                                lower = apply(object, 2, quantile,probs=(1-probs)/2, na.rm=T) ,
                                upper=apply(object, 2, quantile,probs=1-(1-probs)/2, na.rm=T) ,
                                median = apply(object,2,median, na.rm=T)  )

    object_summary = object_summary [order(-object_summary$mean),]
    object_summary$names = rownames(object_summary)
    p = nrow(object_summary)
    object_summary$flag = c(FALSE,object_summary$upper[-1] - object_summary$median[-p]<0)
    object_summary$group = factor(cumsum(object_summary$flag)+1)
  }



  if(  plot.flag ){

    plot_bar = ggplot( object_summary, aes(x = reorder(names, +mean), y = mean) ) +
      geom_bar( stat = "identity" ,fill="grey80")  +
      ylab("Mean(|SHAP|)") + xlab("") +
      geom_errorbar(aes(ymin = lower,ymax= upper,color=group), width = 0.2,linewidth=1)+
      geom_point(aes(y=median),color="blue",size=3) +
      geom_text(aes(label = sprintf("%.3f",round(mean,3))), vjust = -0.3, hjust = -0.3,  size=3.5)+
      theme_bw()  + annotate("text", label = paste0(probs*100,"% quantile") ,x= 1,
                             y = max(object_summary$upper) , vjust = "inward", hjust = "inward") +
      coord_flip(ylim = c(0 , max(object_summary$upper) ))+
      theme(legend.position = "none")

  } else {
    ymin = (max (object_summary$mean) * (-0.1))

    plot_bar = ggplot( object_summary, aes(x = reorder(names, +mean), y = mean) ) +
      geom_bar( stat = "identity" ,fill="grey")  +
      ylab("Mean(|SHAP|)") + xlab("")  +
      geom_text(aes(label = sprintf("%.3f",round(mean,3))), ,  hjust = "inward",  size=3.5)+
      coord_flip(ylim = c(ymin , max(object_summary$mean) ))  + theme_bw()
  }

  if(!is.null(title)) {
    plot_bar = plot_bar + ggtitle(title)
  }
  plot_bar

}


# long format

long_data = function(data, normalize , absolute ) {

  if (isTRUE(absolute)) {
    data_list <- list(raw_data = data.frame(data),
                      mean_data = colMeans(abs(data),na.rm = T)[order(colMeans(abs(data),na.rm = T), decreasing = TRUE)])
  } else if (isFALSE(absolute)) {

    data_list <- list(raw_data = data.frame(data),
                      mean_data = as.data.frame(colMeans(abs(data),na.rm = T )[order(colMeans(abs(data) ,na.rm = T))]))

    names ( data_list $ mean_data ) = "mean"
    data_list $ mean_data $ variable =   rownames (data_list $ mean_data)
    data_list $ mean_data$ rank = seq_len(nrow(data_list $ mean_data))

  }


  names ( data_list $ raw_data) = colnames (data)
  dt  <- data.table::setDT(data_list$raw_data)[, names(data_list$raw_data)[1:dim(data_list$raw_data)[2]], with = FALSE]
  dt [, ID:= .I]


  dt_long <- data.table::melt.data.table(dt, measure.vars = colnames(data_list$raw_data))

  if(isTRUE(normalize)) {
    dt_long  = dt_long  %>% group_by(variable) %>% mutate (normalize = min_max_normalization(value))
  }


  out = list (long = dt_long, mean = data_list $ mean_data)
  return(out)
}


min_max_normalization <- function(x) {
  return ((x - min(x,na.rm = T)) / (max(x,na.rm = T) - min(x,na.rm = T)))
}




# summary plot ------------------

summary_plot = function (object,title,...) {

  plot_summary <- object$long %>% filter(is.na(value)==F) %>%
     ggplot( aes(variable, value)) +
    ggforce::geom_sina(aes(color= normalize*12),
                       method = "counts", maxwidth = 0.7  ) +
    scale_x_discrete(limits = rev(names (object$mean)),
                     labels = rev(names (object$mean))) +
    scale_color_gradient(low="blue", high= "red", breaks=c(0,12), labels=c("Low", "High"),
                         guide = guide_colorbar(barwidth = 0.3, barheight = 12))+
    theme_bw(base_size = 12) + coord_flip()  +
    theme(axis.line.y = element_blank(),axis.ticks.y = element_blank(),
          legend.position="right",legend.title = element_blank() ) +
    labs(y = "SHAP (Impact on model output)", x = "" )

  if(!is.null(title)) {
    plot_summary  =  plot_summary  + ggtitle(title)
  }

  plot_summary

}


# one-hot encoding -> decoded
label_data <- function(data) {

  data = as.data.frame(data)

  var_tmp = uniq_len = factor_value = NULL
  for (i in  names(data) ) {
    temp = unlist  (strsplit(i, split =  "_" ) )
    var_tmp =  c( var_tmp, temp[1] )
    uniq_len = c( uniq_len, length (unique( data [,i])))
    factor_value = c(factor_value,temp[length(temp)] )
  }

  factor_check  =  data.frame(varname =  names(data),
                              var_tmp = var_tmp, uniq_len=uniq_len,
                              factor_value = factor_value)
  factor_check  =  factor_check %>% group_by(var_tmp) %>% mutate(n=n())

  fac = unique(factor_check$var_tmp [factor_check $ uniq_len ==2 & factor_check $ n >=2])

  decoded_data = data
  for (i in fac){
    decoded_data =  decoded_data  %>%  pivot_longer(col=starts_with(i), names_to=i, names_prefix= paste0(i,"."))  %>%
      filter(value==1) %>% select (-value)
  }

  decoded_data = as.data.frame(decoded_data)

  out =list (decoded_data = decoded_data, factor_check = factor_check)
  return (out)
}

