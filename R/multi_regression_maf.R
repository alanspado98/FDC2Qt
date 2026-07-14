#' Multiregression Model for Mean Annual Flow Prediction
#'
#' Computes a linear regression model between mean annual flows and catchment
#' descriptors.
#'
#' @param mean_annual_flows A data frame with two columns: station name
#'   (`character`) and discharge value (`numeric`).
#' @param basin_descriptors A data frame with the same rows as
#'   `mean_annual_flows`, containing multiple columns: station name in the first
#'   (`character`) and morpho‑climatic descriptors in the remaining columns
#'   (`numeric`).
#' @param target A data frame with one row and the same columns as
#'   `basin_descriptors`. Alternatively, it can be a character string containing
#'   the name of a station from `basin_descriptors`. In this case, the station is
#'   considered ungauged and is not included in the model construction.
#' @param m Number of predictors to retain in the model (`numeric`). If not set,
#'   the user can manually define the model by progressively removing
#'   descriptors.
#'
#' @return A list produced by the `lmrob` function, computed for the selected
#'   multiregression model.
#'
#' @examples
#' # Assuming you have mean_annual_flows, basin_descriptors and target
#' # loaded:
#' multi_regression_maf(mean_annual_flows,basin_descriptors,
#'                      target,m = NULL)
#'
#' # This evaluates the multiregression model recursively, allowing the user
#' # to remove descriptors by entering the appropriate column name
#' (`basin_descriptors`).
#'
#Prompts for user removal descriptors
ask_which_descriptor <- function(prompt = "Which descriptor/s?:")
{
  repeat
  {
    response <- trimws(readline(prompt))
    response <- strsplit(response,",")[[1]]
    if(length(which(response %in% names(lm_step$coefficients)[-1]))>0) return(response)
    else cat("Invalid input. Please enter the proper descriptor/s name/s.\n")
  }
}
ask_yes_no <- function(prompt = "Do you want to delete one of these descriptors?[Y/N]:")
{
  repeat
  {
    response <- toupper(trimws(readline(prompt)))
    if(response %in% c("Y", "N")) return(response)
    else cat("Invalid input. Please enter 'Y' or 'N'.\n")
  }
}
#' @export
multi_regression_maf<-function(mean_annual_flows,basin_descriptors,target,m=NULL)
{
  #Common and sorted sections between the 2 datasets
  mean_annual_flows<-mean_annual_flows[mean_annual_flows[,1]%in%basin_descriptors[,1],]
  basin_descriptors<-basin_descriptors[basin_descriptors[,1]%in%mean_annual_flows[,1],]
  basin_descriptors<-basin_descriptors[order(factor(basin_descriptors[,1],levels=mean_annual_flows[,1])),]

  #Save descriptors of the target site
  target_descriptors<-basin_descriptors[basin_descriptors[,1]%in%target,]
  #Removal from MAFs data frame and from descriptors data frame
  if(dim(mean_annual_flows)[1]==dim(basin_descriptors)[1]) #only in case is present
  {
    #(gauged site that we are supposing ungauged)
  mean_annual_flows<-mean_annual_flows[-which(mean_annual_flows[,1]%in%target),]
  basin_descriptors<-basin_descriptors[-which(basin_descriptors[,1]%in%target),]
  }
  #Linear model definition: since in literature a power relationship between MAFs
  # and catchment descriptors has been proven to exist, logarithmic transformation is needed
  #Removing descriptors with at least one negative or NA value
  index_c<-which(sapply(basin_descriptors,function(column) any(column<=0 | is.na(column))))
  if(length(index_c)>0) basin_descriptors <- basin_descriptors[,-index_c]

  user_input_1="Y"
  while(user_input_1=="Y")
  {
    #log transformation
    y<-log(mean_annual_flows[,2])
    x<-log(basin_descriptors[,-1])

    #Best predicting model (according to AIC metric)
    lm_0 <- lm(y ~ 1, data = x)#no predictors, intercept only
    lm_all <- lm(y ~ ., data = x)# all predictors
    lm_step <- step(lm_0,#starting linear model
                    scope = list(lower = lm_0, upper = lm_all),#range of models examined
                    direction="both",#AIC evaluated in removing (backward) or adding
                                     #(forward) each descriptor
                    trace=0)#no information printed during each step

    if(is.null(m))
    {
      print(lm_step$anova)#table for comparing the candidate models
      user_input_1<-ask_yes_no() #to user: descriptor/s removal?
    }else user_input_1<-"N"

    if(user_input_1=="Y")
    {
      user_input_2<-ask_which_descriptor()#to user: which descriptor/s?
      basin_descriptors<-basin_descriptors[,-which(colnames(basin_descriptors)%in%user_input_2)]
      #Descriptor removal
    }
  }
  if(is.null(m)) m<-length(lm_step$coefficients)-1

  formula_user<-as.formula(paste("y ~", paste(names(lm_step$coefficients)[2:(m+1)],collapse = "+")))#symbolic description of the model
  lm_user <- robustbase::lmrob(formula_user,data=x)#linear regression (automatic outlier detection)
  target_descriptors<-target_descriptors[,colnames(target_descriptors)%in%names(lm_step$coefficients)[2:(m+1)]]#target descriptors selected
  target_descriptors<-target_descriptors[order(factor(colnames(target_descriptors),levels=names(lm_step$coefficients)[2:(m+1)]))]
  target_descriptors<-log(target_descriptors)#log transformation
  lm_user$predicted.value<-exp(predict(lm_user,target_descriptors))#append mean annual flow prediction

  return(lm_user)
}
