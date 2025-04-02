#' plot
#'
#' This function uses XXXXXXX
#'
#' @param object explainBART object
#' @param ... add
#' @return Return plot:
#' \item{out}{A \code{plot} plot of Mean (|SHAP|)}
#' @export
plot.explain <-  function(object, average = NULL, type = NULL,  num_post= NULL,
                          plot.flag=TRUE, adjust= FALSE,  probs=0.95,title=NULL,...){

    feature_names = names(object $newdata)
    factor_names = object$factor_names

    length_cate = NULL
    for (i in  factor_names ){
      if(length(which(str_detect(feature_names,i))) > 1  ){
        length_cate = c(  length_cate , feature_names[which(str_detect(feature_names,i))] )
      }
    }

    if(isTRUE(adjust)) {
      # Compute mean and variance of posterior samples
      acomb <- function(...) abind::abind(..., along = 3)
      phis.stats <- foreach(i = feature_names, .combine = "acomb") %do% {
        cbind(rowMeans( object$phis[[i]] ), apply( object$phis[[i]] , MARGIN = 1, FUN = var))
      }
      for (i in seq_len(dim(phis.stats)[1L])) {  # loop through each observation
        err <- object$ fx - sum(phis.stats[i, 1L, ]) - object$fnull
        if( max(phis.stats[i, 2L, ])==0) {  v <- 0  } else {
          v <- (phis.stats[i, 2L, ] / max(phis.stats[i, 2L, ])) * 1e6 }
        adj <- err[i] * (v - (v * sum(v)) / (1 + sum(v)))
        phis.stats[i, 1L, ] <- phis.stats[i, 1L, ] + adj  # adjust Shapley estimate
      }
      phis_adj <- phis.stats[, 1L, ]

      phis_adj = as.data.frame(phis_adj)
      names ( phis_adj) = feature_names


      for ( i in  length_cate ){
        temp_i = which (  object$newdata [,i] !=1)
        phis_adj[temp_i,i] = NA
      }

      local_obs = long_data ( phis_adj, normalize = TRUE, absolute = T)

      if (isTRUE(type =="bar")| is.null (type)){

        Local_mean = as.data.frame(abs( phis_adj ))
        # bar plot
        sample_summary =  local_obs $  long %>% filter (is.na (value)==F  )%>%
          select(c("variable", "value"))%>% group_by(variable) %>% mutate (n=n() ) %>%
          select (c("variable","n"))%>%distinct(.keep_all = T)

        sample_zero = local_obs $  long %>% filter (is.na (value)==F  )%>%
          select(c("variable", "value"))%>% group_by(variable) %>% mutate (n=n() )%>%
          filter(value==0) %>% mutate (zero =n()) %>% distinct(.keep_all = T) %>% select(c("variable", "zero"))

        sample_summary = left_join(sample_summary, sample_zero, by = "variable" )
        sample_summary [which( is.na(sample_summary$zero)), "zero"] = 0

        sample_summary$ percent = paste0(sprintf("%.1f", round(sample_summary$zero/sample_summary$n*100, 1)) ,"%")

        out = bar_plot (Local_mean,probs = NULL, plot.flag = FALSE , title  )

        ymin = - diff(range(out$coordinates$limits$y)) * 0.02

        out = out +
          annotate("text",x= sample_summary$variable, y =  ymin,
                   hjust = 1, label = sample_summary$percent , size = 3.5)


      } else if (isTRUE( type == "bees")){
        #  summary plot
        out = summary_plot (local_obs,title )
        print (out)
      }


    }
    if (is.null(average)){

      stop("If you want to calculate the average based on observations, enter 'obs'.",
           "If you want to calculate the average based on the posterior sample, enter 'post'.",
           "To compare the results based on observation and posterior sample, enter 'both'.",call. = FALSE)

    } else if ( isTRUE(average =="obs")  ){

      for( i in length_cate) {
        temp_i = which(object$newdata[,i] !=1)
        object $ phis [[i]][temp_i,] = NA
      }


      # Average based on observation= #post by #variable
      obs_mean_shap = foreach(i = 1 : length(object $ phis), .combine = 'cbind' ) %do% {
        apply( object $ phis[[i]], 2, function(x) {mean(abs(x),na.rm = T)})
      }
      colnames(obs_mean_shap) = names(object $ phis)


      out = bar_plot (obs_mean_shap, probs, plot.flag= T,title )

      if ( isTRUE(type =="bees") ) {
        message("In Average based on observation,\n a barplot is also provided by entering type='bees'.")
      }

      print (out)


    } else if ( isTRUE( average =="post") ){

      if (is.null(num_post)){


        # Average based on posterior sample = #obs by #variable
        post_mean_shap = foreach(i = 1: length(object $ phis), .combine = 'cbind' ) %do% {
          apply(object $ phis[[i]], 1, function(x) {mean( x) })
        }
        colnames(post_mean_shap) = names(object $phis)

        for ( i in length_cate ){
          temp_i = which ( object$newdata[,i] !=1)
          post_mean_shap[temp_i,i] = NA
        }

        # Change data format & mean values => local_total is List
        local_total = long_data (post_mean_shap, normalize = TRUE, absolute = TRUE)

        if (isTRUE(type =="bar")| is.null (type)){
          # bar plot
          Local_mean = as.data.frame(abs(post_mean_shap))

          sample_summary = local_total $  long %>% filter (is.na (value)==F  )%>%
            select(c("variable", "value"))%>% group_by(variable) %>% mutate (n=n() ) %>%
            select (c("variable","n"))%>%distinct(.keep_all = T)

          sample_zero =  local_total $  long %>% filter (is.na (value)==F  )%>%
            select(c("variable", "value"))%>% group_by(variable) %>% mutate (n=n() )%>%
            filter(value==0) %>% mutate (zero =n()) %>% distinct(.keep_all = T) %>% select(c("variable", "zero"))

          sample_summary = left_join(sample_summary, sample_zero, by = "variable" )
          sample_summary [which( is.na(sample_summary$zero)), "zero"] = 0


          sample_summary$ percent = paste0(sprintf("%.1f", round(sample_summary$zero/sample_summary$n*100, 1)) ,"%")

          out = bar_plot (Local_mean, probs = NULL, plot.flag = FALSE,title  )
          ymin = - diff(range(out$coordinates$limits$y)) * 0.02

          out = out +
            annotate("text",x= sample_summary$variable, y =  ymin,
                     hjust = 1, label = sample_summary$percent , size = 3.5)

          print (out)

        } else if (isTRUE( type == "bees")){
          #  summary plot
          out = summary_plot (local_total,title )
          print (out)

        }
      }


      else if (isFALSE (is.null(num_post)) ){

        # num_post th posterior sample = #obs by #variable
        obs = foreach(i = 1: length(object $ phis), .combine = 'cbind' ) %do% {
          object $ phis[[i]] [,num_post]
        }
        colnames(obs) = names(object $ phis)

        for ( i in length_cate ){
          temp_i = which (  object$newdata [,i] !=1)
          obs[temp_i,i] = NA
        }


        local_obs = long_data (obs, normalize = TRUE, absolute = TRUE)

        if (isTRUE(type == "bees") | is.null (type)){
          # Local SHAP summary plot - i th posterior
          out = summary_plot (local_obs ,title)
          ggpubr::annotate_figure(out, top =  ggpubr::text_grob(paste0("Sample number = ", num_post),
                                                                hjust = 1, x = 0.95,vjust = 1.05, size = 10 ))
        } else if (isTRUE(type == "bar")) {

          out = bar_plot (abs(obs), probs =NULL, plot.flag = FALSE,title )

          sample_summary =  local_obs $  long %>% filter (is.na (value)==F  )%>%
            select(c("variable", "value"))%>% group_by(variable) %>% mutate (n=n() ) %>%
            select (c("variable","n"))%>%distinct(.keep_all = T)

          sample_zero =  local_obs $  long  %>% filter (is.na (value)==F  )%>%
            select(c("variable", "value"))%>% group_by(variable) %>% mutate (n=n() )%>%
            filter(value==0) %>% mutate (zero =n()) %>% distinct(.keep_all = T) %>% select(c("variable", "zero"))

          sample_summary = left_join(sample_summary, sample_zero, by = "variable" )
          sample_summary [which( is.na(sample_summary$zero)), "zero"] = 0

          sample_summary$ percent =  paste0(sprintf("%.1f", round( sample_summary$zero/  sample_summary$n*100, 1)) ,"%")

          ymin = - diff(range(out$coordinates$limits$y)) * 0.02


          ggpubr::annotate_figure(out +  annotate("text",x= sample_summary$variable, y =  ymin,
                                                  hjust = 1, label = sample_summary$percent , size = 3.5),
                                  top =  ggpubr::text_grob(paste0("Sample number = ", num_post),
                                                           hjust = 1, x = 0.95,vjust = 1.05, size = 10 ))

        }
      }



    } else if( isTRUE(average=="both") ) {

      for( i in length_cate) {
        temp_i = which(object$newdata[,i] !=1)
        object $ phis[[i]][temp_i,] = NA
      }

      # Average based on observation= #post by #variable
      obs_mean_shap = as.data.frame(foreach(i = 1 : length(object $ phis), .combine = 'cbind' ) %do% {
        apply(object $ phis[[i]], 2, function(x) {mean(abs(x), na.rm= T)})
      })
      colnames(obs_mean_shap) = names(object $ phis)


      mean_obs = colMeans(obs_mean_shap) [order( colMeans(obs_mean_shap), decreasing = TRUE) ]
      mean_obs = as.data.frame(mean_obs)
      names (mean_obs) = "mean"
      mean_obs$ Average ="Observation"
      mean_obs$ variable= obs_variable = rownames (mean_obs)

      obs_mean_shap $ Average ="Observation"
      obs_mean_shap $ ID = rownames(obs_mean_shap)



      # Average based on posterior sample = #obs by #variable

      post_mean_shap = as.data.frame(foreach(i = 1: length(object $ phis), .combine = 'cbind' ) %do% {
        apply(object $ phis[[i]], 1, function(x) {abs(mean( x ) )})
      })
      colnames(post_mean_shap) = names(object $ phis)

      for ( i in length_cate ){
        temp_i = which ( object$newdata[,i] !=1)
        post_mean_shap[temp_i,i] = NA
      }

      mean_post = colMeans(post_mean_shap,na.rm = T) [order( colMeans(post_mean_shap,na.rm = T ), decreasing = TRUE) ]
      mean_post = as.data.frame(mean_post)
      names (mean_post) = "mean"
      mean_post $ Average = "Samples"
      mean_post$ variable = rownames (mean_post)

      post_mean_shap$ Average ="Samples"
      post_mean_shap $ ID = rownames(post_mean_shap)


      total_mean_shap = rbind (obs_mean_shap, post_mean_shap )

      total_long = reshape2::melt(total_mean_shap, id =c('ID',"Average"))
      total_long $ Average_re = forcats::fct_rev( total_long $ Average )

      mean_data =  rbind (mean_obs, mean_post)
      mean_data  $ Average_re = forcats::fct_rev( mean_data  $ Average)

      data_list <- list(total_long = total_long, mean_data = mean_data  )


      if ( is.null (type) | isTRUE(type == "bees")) {

        ymax = max (data_list $ total_long$value, na.rm = T) * (1+0.001)
        ymin = - diff(range(data_list $ total_long$value))*0.05

        left_p =  data_list $ total_long%>% filter (is.na (value)==F & Average_re=="Observation")%>%
          ggplot(  aes(x= variable, y=value, fill= Average_re )) +
          geom_boxplot() +coord_flip() + theme_bw(base_size = 11) +
          labs(fill = "", x="",y="Mean(|SHAP|)") +
          scale_x_discrete(limits = rev(obs_variable),
                           labels = rev(obs_variable))+
          scale_fill_manual(values=c(  "#F8766D"))+
          theme(legend.position = "bottom") + ylim(ymin, ymax )

        sample_summary = data_list $ total_long %>% filter (is.na (value)==F & Average_re=="Samples")%>%
          select(c("variable", "value"))%>% group_by(variable) %>% mutate (n=n() ) %>%
          select (c("variable","n"))%>%distinct(.keep_all = T)

        sample_zero = data_list $ total_long %>% filter (is.na (value)==F & Average_re=="Samples" )%>%
          select(c("variable", "value"))%>% group_by(variable) %>% mutate (n=n() )%>%
          filter(value==0) %>% mutate (zero =n()) %>% distinct(.keep_all = T) %>% select(c("variable", "zero"))

        sample_summary = left_join(sample_summary, sample_zero, by = "variable" )
        sample_summary [which( is.na(sample_summary$zero)), "zero"] = 0

        sample_summary$ percent = paste0(sprintf("%.1f", round(sample_summary$zero/sample_summary$n*100, 1))  ,"%")

        right_p =  data_list $ total_long%>% filter (is.na (value)==F & Average_re=="Samples")%>%
          ggplot(  aes(x= variable, y=value, fill= Average_re )) +
          geom_boxplot() +coord_flip() + theme_bw(base_size = 11) +
          labs(fill = "", x="",y="Mean(|SHAP|)") +
          scale_x_discrete(limits = rev(obs_variable),
                           labels = rev(obs_variable))+
          scale_fill_manual(values=c(  "#00BFC4"))+
          theme(legend.position = "bottom") + ylim(ymin, ymax )

        right_p = right_p + annotate("text",x= sample_summary$variable,
                                     y = - diff(range(data_list $ total_long$value))*0.02 ,
                                     hjust = 1, label = sample_summary$percent , size = 2.5)


        gridExtra ::grid.arrange(left_p, right_p, ncol= 1, top=title  )

        if ( isTRUE(type =="bees") ) {
          message("In Average = 'both',\n a baxplot is also provided by entering type='bees'.")
        }

      } else if (isTRUE(type == "bar")  ){
        out = ggplot( data_list $ mean_data, aes(x= variable, y= mean, fill =  Average_re )) +
          geom_bar( stat = "identity" , position="dodge") +coord_flip() + theme_bw(base_size = 11) +
          scale_fill_manual(values=c( "#00BFC4","#F8766D"))+
          labs(fill = "Criteria for average ", x="",y="Mean(|SHAP|)") +
          scale_x_discrete(limits = rev(obs_variable), labels = rev(obs_variable)) +
          ylim(c(0,max( data_list $ mean_data$mean)  ))+
          geom_text(aes( label = sprintf("%.3f",round(mean,3)) ),
                    hjust ="inward",  size= 3.5, position = position_dodge(width=0.9))+
          theme(legend.position = "bottom")

        if(!is.null(title)) {
          out = out + ggtitle(title)
        }

        print(out)

      }

    }

  }




