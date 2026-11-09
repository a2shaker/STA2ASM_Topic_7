# STA2ASM Week 7 - Computational inference: Bootstrap Explorer
# Requires shiny and rquest. Uses Base R graphics and the established Week 4 palette.
library(shiny)
library(rquest)
if (!requireNamespace("rquest", quietly = TRUE)) {
  stop("Please install the rquest package: install.packages(\"rquest\")")
}

PAL <- list(teal="#236A6E", blue="#174F8A", rose="#C65D7B", purple="#7030A0",
            original="#E9CEDF", bootstrap="#CDE4E3", green="#548235",
            coral="#C0504D", pale="#eef5f5")
fmt <- function(x, digits=3) ifelse(is.finite(x), formatC(x,format="f",digits=digits), "Not defined")
math_label <- function(label, symbol) withMathJax(HTML(paste0(label, " \\(", symbol, "\\)")))
# Inline mathematical notation in dynamically generated Shiny text.
mj <- function(value) withMathJax(HTML(value))
choices <- c("Normal"="norm", "Log-normal"="lnorm", "Exponential"="exp",
             "Uniform"="unif", "Student's t"="t", "Chi-squared"="chisq",
             "Bernoulli"="bern", "Binomial"="bin", "Poisson"="pois",
             "Negative binomial (total trials)"="nbt",
             "Negative binomial (mean-dispersion)"="nba")
stats <- c("Mean"="mean", "Median"="median", "Selected quantile"="quantile",
           "Standard deviation"="sd", "Variance"="variance", "Coefficient of variation (CV)"="cv")

specification <- function(dist, mu=0, sigma=1, df=5, p=.3, trials=20,
                          lambda=2.6, size=3, nbmu=7, theta=3, lo=0, hi=1) {
  switch(dist,
         norm=list(label="Normal",gen=function(n) rnorm(n,mu,sigma),q=function(p) qnorm(p,mu,sigma),mean=mu,var=sigma^2),
         lnorm=list(label="Log-normal",gen=function(n) rlnorm(n,mu,sigma),q=function(p) qlnorm(p,mu,sigma),mean=exp(mu+sigma^2/2),var=(exp(sigma^2)-1)*exp(2*mu+sigma^2)),
         exp=list(label="Exponential",gen=function(n) rexp(n,lambda),q=function(p) qexp(p,lambda),mean=1/lambda,var=1/lambda^2),
         unif=list(label="Uniform",gen=function(n) runif(n,lo,hi),q=function(p) qunif(p,lo,hi),mean=(lo+hi)/2,var=(hi-lo)^2/12),
         t=list(label="Student's t",gen=function(n) rt(n,df),q=function(p) qt(p,df),mean=if(df>1) 0 else NA_real_,var=if(df>2) df/(df-2) else Inf),
         chisq=list(label="Chi-squared",gen=function(n) rchisq(n,df),q=function(p) qchisq(p,df),mean=df,var=2*df),
         bern=list(label="Bernoulli",gen=function(n) rbinom(n,1,p),q=function(q) qbinom(q,1,p),mean=p,var=p*(1-p)),
         bin=list(label="Binomial",gen=function(n) rbinom(n,trials,p),q=function(q) qbinom(q,trials,p),mean=trials*p,var=trials*p*(1-p)),
         pois=list(label="Poisson",gen=function(n) rpois(n,lambda),q=function(q) qpois(q,lambda),mean=lambda,var=lambda),
         nbt=list(label="Negative binomial (total trials)",gen=function(n) rnbinom(n,size,p)+size,q=function(q) qnbinom(q,size,p)+size,mean=size/p,var=size*(1-p)/p^2),
         nba=list(label="Negative binomial (mean-dispersion)",gen=function(n) rnbinom(n,size=theta,mu=nbmu),q=function(q) qnbinom(q,size=theta,mu=nbmu),mean=nbmu,var=nbmu+nbmu^2/theta)
  )
}
statistic <- function(x, what, prob=.75) {
  switch(what,mean=mean(x),median=median(x),quantile=unname(quantile(x,prob)),
         sd=sd(x),variance=var(x),cv=if(abs(mean(x))<1e-10) NA_real_ else 100*sd(x)/mean(x))
}
truth <- function(s, what, prob=.75) {
  switch(what,mean=s$mean,median=s$q(.5),quantile=s$q(prob),
         sd=sqrt(s$var),variance=s$var,
         cv=if(!is.finite(s$mean)||abs(s$mean)<1e-10||!is.finite(s$var)) NA_real_ else 100*sqrt(s$var)/s$mean)
}
boot_stats <- function(x, what, B, prob=.75) {
  n <- length(x); out <- numeric(B)
  for(b in seq_len(B)) out[b] <- statistic(sample(x,n,replace=TRUE),what,prob)
  out
}
cox_interval <- function(x, level=.95) {
  if(length(x)<2 || any(x<=0)) return(c(NA_real_,NA_real_))
  y<-log(x); n<-length(y); v<-var(y)
  center<-mean(y)+v/2; se<-sqrt(v/n+v^2/(2*(n-1)))
  exp(center+c(-1,1)*qt(1-(1-level)/2,n-1)*se)
}
# rquest asymptotic quantile CI; type 7 matches Base R quantile().
# Discrete or very small samples can yield undefined or non-finite results.
rquest_interval <- function(x, what, level, prob=.75) {
  if (!(what %in% c("median", "quantile")) || length(x) < 3 ||
      any(!is.finite(x))) return(c(NA_real_, NA_real_))
  ans <- tryCatch(suppressWarnings(
    if (what == "median")
      rquest::q.test(x, measure="median", quantile.type=7,
                     conf.level=level)$conf.int
    else
      rquest::q.test(x, u=prob, quantile.type=7,
                     conf.level=level)$conf.int
  ), error=function(e) c(NA_real_,NA_real_))
  ans <- as.numeric(ans)
  if (length(ans)!=2 || any(!is.finite(ans)) || ans[1]>ans[2])
    return(c(NA_real_,NA_real_))
  ans
}

ci_methods <- function(x, b, what, level, dist, prob=.75) {
  alpha<-1-level
  result<-list("Bootstrap percentile"=unname(quantile(b,c(alpha/2,1-alpha/2),na.rm=TRUE)))
  if (what %in% c("median", "quantile"))
    result[["rquest quantile interval"]] <- rquest_interval(x,what,level,prob)
  if(what=="mean" && length(x)>1 && all(is.finite(x))) {
    # A constant sample has a zero-width t interval; t.test() otherwise errors.
    if(sd(x)==0) result[["Theoretical t-interval"]] <- rep(mean(x),2)
    else result[["Theoretical t-interval"]] <- unname(t.test(x,conf.level=level)$conf.int[1:2])
  }
  if(what=="median" && dist=="lnorm" && length(x)>1 && all(is.finite(x)) && all(x>0)) {
    y <- log(x)
    margin <- qt(1-alpha/2,df=length(x)-1)*sd(y)/sqrt(length(x))
    result[["Back-transformed theoretical t-interval"]] <- exp(mean(y)+c(-1,1)*margin)
  }
  if(what=="mean" && dist %in% c("lnorm","exp","chisq"))
    result[[if(dist=="lnorm") "Modified Cox" else "Modified Cox (assumed log-normal)"]] <- cox_interval(x,level)
  result
}

parameter_ui <- function(prefix, dist) {
  id<-function(x) paste0(prefix,"_",x)
  switch(dist,
         norm=tagList(numericInput(id("mu"),math_label("Mean", "\\mu"),0),numericInput(id("sigma"),math_label("Standard deviation", "\\sigma"),1,min=.1,max=10,step=.1)),
         lnorm=tagList(numericInput(id("mu"),math_label("Mean on log scale", "\\mu"),0),numericInput(id("sigma"),math_label("SD on log scale", "\\sigma"),1,min=.1,max=3,step=.1)),
         exp=numericInput(id("lambda"),math_label("Rate", "\\lambda"),1,min=.05,max=20,step=.1),
         unif=tagList(numericInput(id("lo"),"Lower bound",0),numericInput(id("hi"),"Upper bound",10)),
         t=numericInput(id("df"),"Degrees of freedom",5,min=1,max=100,step=1),
         chisq=numericInput(id("df"),"Degrees of freedom",5,min=1,max=100,step=1),
         bern=sliderInput(id("p"),"Probability of success",.01,.99,.3,step=.01),
         bin=tagList(numericInput(id("trials"),"Number of trials",20,min=1,max=300,step=1),sliderInput(id("p"),"Probability of success",.01,.99,.3,step=.01)),
         pois=numericInput(id("lambda"),math_label("Expected count", "\\lambda"),3,min=.05,max=100,step=.5),
         nbt=tagList(numericInput(id("size"),"Required successes",3,min=1,max=100,step=1),sliderInput(id("p"),"Probability of success",.01,.99,.3,step=.01)),
         nba=tagList(numericInput(id("nbmu"),"Mean",7,min=.05,max=100,step=.5),numericInput(id("theta"),math_label("Dispersion", "\\theta"),3,min=.1,max=100,step=.5))
  )
}

sample_controls <- function(prefix, include_B=TRUE) {
  tagList(selectInput(paste0(prefix,"_dist"),"Population distribution",choices,selected="lnorm"),
          uiOutput(paste0(prefix,"_params")),
          selectInput(paste0(prefix,"_stat"),"Statistic",stats,selected="mean"),
          conditionalPanel(sprintf("input.%s_stat == 'quantile'",prefix),
                           sliderInput(paste0(prefix,"_q"),"Quantile probability",.05,.95,.75,step=.05)),
          numericInput(paste0(prefix,"_n"),math_label("Original sample size", "n"),30,min=5,max=500,step=5),
          if(include_B) selectInput(paste0(prefix,"_B"),math_label("Bootstrap resamples", "B"),c("100"=100,"500"=500,"1,000"=1000,"2,000"=2000,"5,000"=5000,"10,000"=10000),selected=2000),
          sliderInput(paste0(prefix,"_level"),"Confidence level",.8,.99,.95,step=.01,ticks=FALSE))
}

ui <- fluidPage(
  tags$head(tags$style(HTML('
    body {font-family:"Helvetica Neue",Arial,sans-serif;color:#333;background:white;}
    h1,h2,h3,h4 {color:#236A6E;font-weight:500;}
    .well {background:#eef5f5;border:1px solid #cfdfe0;border-left:4px solid #4d8386;}
    .navbar-default {background:white;border:0;border-bottom:1px solid #dce6e6;border-radius:0;box-shadow:none;margin-bottom:16px;}
    .navbar-default .navbar-nav>li>a {color:#4a5b5c;padding:12px 14px 9px;border-bottom:3px solid transparent;}
    .navbar-default .navbar-nav>.active>a,.navbar-default .navbar-nav>.active>a:hover {color:#236A6E;font-weight:600;background:transparent;border-bottom:3px solid #236A6E;}
    .btn-primary {background:#236A6E;border-color:#236A6E;}
    .result {background:#eef5f5;border-left:4px solid #236A6E;padding:12px 15px;margin:12px 0;}
    .note {background:#F3F5F7;border-left:3px solid #AAB8C2;padding:10px 14px;margin-top:12px;color:#4C5860;}
    .subtle {color:#58696a;font-size:13px;}
    pre {white-space:pre-wrap;}
  '))),
  titlePanel("STA2ASM - Week 7: Computational Inference"),
  withMathJax(),
  p(class="subtle","Explore bootstrap resampling, confidence intervals, coverage and null-centred hypothesis testing."),
  navbarPage(title=NULL,id="main_tab",collapsible=TRUE,
             tabPanel("Explore the Bootstrap",value="explore",
                      fluidRow(column(3,wellPanel(sample_controls("ex"),
                                                  actionButton("ex_new","Generate new sample",class="btn-primary"),
                                                  actionButton("ex_resample","Resample fixed data"),
                                                  actionButton("ex_reset","Restore reading defaults"))),
                               column(9,h2("Explore the Bootstrap"),
                                      p("The original sample stays fixed when you resample. Generate a new sample to investigate sampling variability."),
                                      fluidRow(column(6,plotOutput("ex_original",height="365px")),
                                               column(6,plotOutput("ex_bootstrap",height="365px"))),
                                      uiOutput("ex_results"),tags$hr(),
                                      h3("Example individual bootstrap samples"),
                                      uiOutput("ex_draws"),
                                      NULL))),
             tabPanel("Compare Confidence Intervals",value="compare",
                      fluidRow(column(3,wellPanel(sample_controls("co"),
                                                  actionButton("co_new","Generate new sample",class="btn-primary"),
                                                  actionButton("co_resample","Resample fixed data"),
                                                  actionButton("co_reset","Restore reading defaults"))),
                               column(9,h2("Compare Confidence Intervals"),
                                      mj("<p>Compare bootstrap percentile intervals with the theoretical \\(t\\)-interval for the mean and the modified Cox interval where applicable.</p>"),
                                      plotOutput("co_plot",height="410px"),uiOutput("co_results"),
                                      conditionalPanel("input.co_stat == 'median' || input.co_stat == 'quantile'",
                                                       div(class="subtle",tags$code("rquest")," intervals use quantile type 7. Results for discrete populations may be unreliable or unavailable.")),
                                      conditionalPanel("input.co_stat != 'mean' && !(input.co_stat == 'median' && input.co_dist == 'lnorm')",
                                                       div(class="subtle",mj("Theoretical \\(t\\)-intervals for the mean and modified Cox intervals are available when Mean is selected."))),
                                      conditionalPanel("input.co_stat == 'mean' && (input.co_dist == 'exp' || input.co_dist == 'chisq')",
                                                       div(class="subtle","Modified Cox is shown as an exploratory assumed-log-normal comparison; exponential and chi-squared populations are not log-normal."))))),
             tabPanel("Investigate Coverage",value="coverage",
                      fluidRow(column(3,wellPanel(sample_controls("cv"),
                                                  selectInput("cv_m","Independent simulated samples",c("100"=100,"250"=250,"500"=500,"1,000"=1000),selected=250),
                                                  selectInput("cv_show","Intervals to display",c("20"=20,"50"=50,"100"=100),selected=50),
                                                  actionButton("cv_run","Run coverage simulation",class="btn-primary"),
                                                  actionButton("cv_reset","Restore reading defaults"))),
                               column(9,h2("Investigate Coverage"),
                                      p("Generate independent samples from a known population, calculate confidence intervals and check how often they contain the true parameter."),
                                      uiOutput("cv_results"),uiOutput("cv_method_ui"),plotOutput("cv_plot",height="420px"),
                                      div(class="subtle","Teal intervals contain the true parameter; coral intervals miss it. Coverage is based on all simulated samples, not just those displayed."),
                                      NULL))),
             tabPanel("Bootstrap Hypothesis Testing",value="test",
                      fluidRow(column(3,wellPanel(
                        selectInput("ht_dist","Population distribution",choices,selected="norm"),uiOutput("ht_params"),
                        numericInput("ht_n",math_label("Original sample size", "n"),30,min=5,max=500,step=5),
                        selectInput("ht_B",math_label("Null-centred bootstrap resamples", "B"),c("500"=500,"1,000"=1000,"2,000"=2000,"5,000"=5000),selected=2000),
                        numericInput("ht_mu0",math_label("Hypothesised mean", "\\mu_0"),0,step=.5),
                        actionButton("ht_new","Generate new sample",class="btn-primary"),
                        actionButton("ht_resample","Resample under H₀"),
                        actionButton("ht_reset","Restore reading defaults"))),
                        column(9,h2("Bootstrap Hypothesis Testing"),
                               uiOutput("ht_hypotheses"),
                               mj("<p>For a two-sided one-sample mean test, shift observations to the hypothesised mean \\(\\mu_0\\), then compare the observed \\(t\\)-statistic with the null-centred bootstrap distribution.</p>"),
                               plotOutput("ht_plot",height="400px"),uiOutput("ht_results"),
                               div(class="note",mj("Rose bars show bootstrap test statistics at least as extreme as the observed statistic. The theoretical \\(t\\)-test \\(p\\)-value is the area under the blue curve beyond the same thresholds.")))))
  )
)

server <- function(input,output,session) {
  get_spec <- function(prefix,dist) {
    v<-function(name,default) {x<-input[[paste0(prefix,"_",name)]];if(is.null(x)||!is.finite(x)) default else x}
    specification(dist,mu=v("mu",0),sigma=max(.01,v("sigma",1)),df=max(1,round(v("df",5))),
                  p=min(.999,max(.001,v("p",.3))),trials=max(1,round(v("trials",20))),
                  lambda=max(.001,v("lambda",2.6)),size=max(1,round(v("size",3))),
                  nbmu=max(.001,v("nbmu",7)),theta=max(.01,v("theta",3)),
                  lo=v("lo",0),hi=max(v("lo",0)+.001,v("hi",10)))
  }
  for(prefix in c("ex","co","cv","ht")) {
    local({pfx<-prefix
    output[[paste0(pfx,"_params")]]<-renderUI(parameter_ui(pfx,input[[paste0(pfx,"_dist")]]))
    })
  }
  spec_for<-function(prefix) get_spec(prefix,input[[paste0(prefix,"_dist")]])
  q_for<-function(prefix) {q<-input[[paste0(prefix,"_q")]];if(is.null(q)) .75 else q}
  
  # Original sample and bootstrap replicates are separate reactive events.
  ex_sample<-eventReactive(list(input$ex_new,input$ex_dist,input$ex_mu,input$ex_sigma,
                                input$ex_df,input$ex_p,input$ex_trials,input$ex_lambda,
                                input$ex_size,input$ex_nbmu,input$ex_theta,input$ex_lo,input$ex_hi,input$ex_n), {
                                  spec_for("ex")$gen(as.integer(input$ex_n))
                                },ignoreNULL=FALSE)
  ex_boot<-eventReactive(list(input$ex_resample,ex_sample(),input$ex_B,input$ex_stat,input$ex_q),{
    boot_stats(ex_sample(),input$ex_stat,as.integer(input$ex_B),q_for("ex"))
  },ignoreNULL=FALSE)
  stat_name <- function(what, q=.75) {
    switch(what, mean="mean", median="median",
           quantile=sprintf("%.0fth percentile",100*q),
           sd="standard deviation", variance="variance",
           cv="coefficient of variation (CV)")
  }
  output$ex_original<-renderPlot({
    x<-ex_sample(); what<-input$ex_stat; q<-q_for("ex")
    oldmar<-par(mar=c(5,4,7,1)); on.exit(par(oldmar))
    h<-hist(x,plot=FALSE)
    hist(x,breaks=h$breaks,col=PAL$original,border="black",main="",
         xlab="Individual observation",ylab="Frequency",las=1,bty="l")
    if(what %in% c("mean","median","quantile")) {
      v<-statistic(x,what,q)
      abline(v=v,col=PAL$teal,lwd=2.5)
      legend("topright",paste("Sample",stat_name(what,q)),col=PAL$teal,
             lty=1,lwd=2.5,bty="n",cex=.9,
             inset=c(0,-.045),xpd=NA)
    }
    title(main="Original sample",line=5)
  })
  output$ex_bootstrap<-renderPlot({
    oldmar<-par(mar=c(5,4,7,1)); on.exit(par(oldmar))
    x<-ex_boot(); req(all(is.finite(x)))
    h<-hist(x,plot=FALSE)
    ci<-quantile(x,c((1-input$ex_level)/2,1-(1-input$ex_level)/2))
    observed<-statistic(ex_sample(),input$ex_stat,q_for("ex"))
    target<-truth(spec_for("ex"),input$ex_stat,q_for("ex"))
    marks<-c(observed,ci,if(is.finite(target)) target)
    lim<-range(c(h$breaks,marks));pad<-max(diff(lim)*.06,1e-6)
    hist(x,breaks=h$breaks,freq=FALSE,col=PAL$bootstrap,border="black",
         main="",xlab=paste("Bootstrap",stat_name(input$ex_stat,q_for("ex"))),
         ylab="Density",las=1,bty="l",xlim=lim+c(-pad,pad),ylim=c(0,max(h$density)*1.35))
    # Stop reference lines below the legend, which sits in the upper-right.
    ymax<-max(h$density)*1.06
    segments(observed,0,observed,ymax,col=PAL$teal,lwd=2.5)
    segments(ci,0,ci,ymax,col=PAL$rose,lty=2,lwd=2)
    if(is.finite(target)) segments(target,0,target,ymax,col=PAL$blue,lty=3,lwd=2)
    legend("topright",legend=c("Observed statistic",
                               sprintf("%.0f%% percentile CI",100*input$ex_level),
                               if(is.finite(target)) "Population parameter"),
           col=c(PAL$teal,PAL$rose,if(is.finite(target)) PAL$blue),
           lty=c(1,2,if(is.finite(target)) 3),lwd=c(2.5,2,if(is.finite(target)) 2),
           bty="n",cex=.85,inset=c(0,-.045),xpd=NA)
    title(main="Bootstrap distribution",line=5)
  })
  output$ex_results<-renderUI({
    b<-ex_boot(); x<-ex_sample(); what<-input$ex_stat; q<-q_for("ex")
    ci<-quantile(b,c((1-input$ex_level)/2,1-(1-input$ex_level)/2))
    name<-stat_name(what,q)
    suffix<-if(what=="cv") "%" else ""
    value<-function(z) paste0(fmt(z),suffix)
    div(class="result",
        div(style="display:flex;flex-wrap:wrap;gap:30px;margin-bottom:10px;",
            div(strong(paste("Sample",name)),tags$br(),
                tags$span(style="font-size:20px;",value(statistic(x,what,q)))),
            div(strong(paste("True population",name)),tags$br(),
                tags$span(style="font-size:20px;",value(truth(spec_for("ex"),what,q))))),
        strong(sprintf("%.0f%% bootstrap percentile confidence interval",100*input$ex_level)),
        tags$br(),sprintf("(%s, %s)",value(ci[1]),value(ci[2])),tags$br(),
        tags$span(class="subtle",paste("Bootstrap standard error:",value(sd(b)))))
  })
  output$ex_draws<-renderUI({
    x<-ex_sample(); n<-length(x)
    # Fixed illustrative resamples based on the current original data.
    tagList(lapply(seq_len(3),function(b){
      inds<-sample.int(n,n,replace=TRUE)
      div(style="margin-bottom:12px;",strong(paste("Bootstrap sample",b)),tags$br(),
          div(style="display:flex;flex-wrap:wrap;gap:6px;margin-top:6px;",
              lapply(seq_len(min(n,20)),function(i)
                tags$span(style="display:inline-block;min-width:43px;text-align:center;padding:7px 9px;border-radius:6px;background:#F0DCE8;border:1px solid #c9c2d8;font-weight:600;",fmt(x[inds[i]],1)))),
          if(n>20) span(class="subtle",sprintf("Showing first 20 of %d observations",n)))
    }))
  })
  
  co_sample<-eventReactive(list(input$co_new,input$co_dist,input$co_mu,input$co_sigma,
                                input$co_df,input$co_p,input$co_trials,input$co_lambda,
                                input$co_size,input$co_nbmu,input$co_theta,input$co_lo,input$co_hi,input$co_n),{
                                  spec_for("co")$gen(as.integer(input$co_n))
                                },ignoreNULL=FALSE)
  co_boot<-eventReactive(list(input$co_resample,co_sample(),input$co_B,input$co_stat,input$co_q),{
    boot_stats(co_sample(),input$co_stat,as.integer(input$co_B),q_for("co"))
  },ignoreNULL=FALSE)
  co_ci<-reactive(ci_methods(co_sample(),co_boot(),input$co_stat,input$co_level,input$co_dist,q_for("co")))
  output$co_plot<-renderPlot({
    ci<-co_ci(); target<-truth(spec_for("co"),input$co_stat,q_for("co"))
    vals<-unlist(ci); vals<-vals[is.finite(vals)]
    req(length(vals)>0)
    if(is.finite(target)) vals<-c(vals,target)
    xr<-range(vals); pad<-max(diff(xr)*.12,.01)
    k<-length(ci); oldmar<-par(mar=c(5,20,4,2)); on.exit(par(oldmar))
    plot(c(xr[1]-pad,xr[2]+pad),c(.5,k+.5),type="n",yaxt="n",las=1,
         xlab="Population parameter / confidence interval",ylab="",main="Confidence interval comparison",bty="l")
    axis(2,at=seq_len(k),labels=names(ci),las=1,cex.axis=.85)
    cols<-c(PAL$rose,PAL$blue,PAL$purple,PAL$coral)
    for(i in seq_len(k)) if(all(is.finite(ci[[i]]))) {
      segments(ci[[i]][1],i,ci[[i]][2],i,col=cols[(i-1)%%length(cols)+1],lwd=4)
      points(ci[[i]],rep(i,2),pch=16,col=cols[(i-1)%%length(cols)+1],cex=.9)
    }
    if(is.finite(target)) abline(v=target,col=PAL$teal,lty=2,lwd=2)
    if(is.finite(target)) legend("topright","True population parameter",col=PAL$teal,lty=2,lwd=2,bty="n",cex=.8)
  })
  output$co_results<-renderUI({
    ci<-co_ci(); target<-truth(spec_for("co"),input$co_stat,q_for("co"))
    div(class="result",strong(paste("True population parameter:",fmt(target))),tags$br(),
        lapply(names(ci),function(name){
          v<-ci[[name]]
          label<-if(name=="Back-transformed theoretical t-interval") "Back-transformed theoretical \\(t\\)-interval" else if(name=="Theoretical t-interval") "Theoretical \\(t\\)-interval" else if(name=="rquest quantile interval") "<code>rquest</code> quantile interval" else name
          tagList(mj(if(all(is.finite(v))) sprintf("%s: (%s, %s)",label,fmt(v[1]),fmt(v[2])) else paste0(label,": unavailable for this sample")),tags$br())
        }))
  })
  
  cv_data<-eventReactive(input$cv_run,{
    s<-spec_for("cv");n<-as.integer(input$cv_n);B<-as.integer(input$cv_B)
    m<-as.integer(input$cv_m);what<-input$cv_stat;q<-q_for("cv")
    target<-truth(s,what,q);req(is.finite(target))
    method_names<-c("Bootstrap percentile",if(what %in% c("median","quantile")) "rquest quantile interval",if(what=="mean") "Theoretical t-interval",
                    if(what=="median" && input$cv_dist=="lnorm") "Back-transformed theoretical t-interval",
                    if(what=="mean" && input$cv_dist %in% c("lnorm","exp","chisq"))
                      if(input$cv_dist=="lnorm") "Modified Cox" else "Modified Cox (assumed log-normal)")
    out<-lapply(method_names,function(z) matrix(NA_real_,m,2));names(out)<-method_names
    withProgress(message="Simulating coverage",value=0,{
      for(i in seq_len(m)) {
        x<-s$gen(n); b<-boot_stats(x,what,B,q)
        ci<-ci_methods(x,b,what,input$cv_level,input$cv_dist,q)
        for(name in names(out)) if(!is.null(ci[[name]])) out[[name]][i,]<-ci[[name]]
        if(i%%max(1,floor(m/20))==0) incProgress(1/20)
      }
    })
    list(intervals=out,target=target,m=m,level=input$cv_level)
  },ignoreNULL=FALSE)
  output$cv_results<-renderUI({
    d<-cv_data();div(class="result",strong(sprintf("Nominal coverage: %.0f%%; true parameter: %s",100*d$level,fmt(d$target))),tags$br(),
                     lapply(names(d$intervals),function(name){a<-d$intervals[[name]];ok<-is.finite(a[,1])&is.finite(a[,2]);
                     cover<-if(any(ok)) mean(a[ok,1]<=d$target & d$target<=a[ok,2]) else NA_real_
                     width<-if(any(ok)) mean(a[ok,2]-a[ok,1]) else NA_real_
                     label<-if(name=="Theoretical t-interval") "Theoretical \\(t\\)-interval" else if(name=="Back-transformed theoretical t-interval") "Back-transformed theoretical \\(t\\)-interval" else if(name=="rquest quantile interval") "<code>rquest</code> quantile interval" else name
                     tagList(mj(sprintf("%s: coverage %s%%; mean width %s",label,fmt(100*cover,1),fmt(width))),tags$br())}))
  })
  # Select which method to visualise without regenerating the simulated samples.
  output$cv_method_ui <- renderUI({
    d <- cv_data()
    methods <- names(d$intervals)
    current <- isolate(input$cv_method)
    if (is.null(current) || !(current %in% methods)) current <- methods[1]
    selectInput("cv_method", "Confidence interval method to display",
                choices = methods, selected = current)
  })
  output$cv_plot<-renderPlot({
    d<-cv_data()
    method<-input$cv_method
    if(is.null(method) || !(method %in% names(d$intervals))) method<-names(d$intervals)[1]
    a<-d$intervals[[method]]
    k<-min(as.integer(input$cv_show),nrow(a))
    a<-a[seq_len(k),,drop=FALSE]
    ok<-is.finite(a[,1])&is.finite(a[,2]); req(any(ok))
    lim<-range(c(a[ok,],d$target));pad<-max(diff(lim)*.06,.01)
    plot(c(lim[1]-pad,lim[2]+pad),c(.5,k+.5),type="n",las=1,
         xlab="Parameter value",ylab="Simulated sample",
         main=paste(method,"confidence intervals"),bty="l")
    abline(v=d$target,col=PAL$blue,lty=2,lwd=2)
    covered<-ok & a[,1]<=d$target & d$target<=a[,2]
    cols<-ifelse(covered,PAL$teal,PAL$coral)
    segments(a[ok,1],which(ok),a[ok,2],which(ok),col=cols[ok],lwd=2)
    points(rowMeans(a[ok,,drop=FALSE]),which(ok),pch=16,cex=.5,col=cols[ok])
    legend("topright",c("Contains true parameter","Misses true parameter","True parameter"),
           col=c(PAL$teal,PAL$coral,PAL$blue),lty=c(1,1,2),lwd=2,bty="n",cex=.8)
  })
  
  ht_sample<-eventReactive(list(input$ht_new,input$ht_dist,input$ht_mu,input$ht_sigma,
                                input$ht_df,input$ht_p,input$ht_trials,input$ht_lambda,
                                input$ht_size,input$ht_nbmu,input$ht_theta,input$ht_lo,input$ht_hi,input$ht_n),{
                                  spec_for("ht")$gen(as.integer(input$ht_n))
                                },ignoreNULL=FALSE)
  ht_null<-eventReactive(list(input$ht_resample,ht_sample(),input$ht_B,input$ht_mu0),{
    x<-ht_sample();n<-length(x);mu0<-input$ht_mu0;B<-as.integer(input$ht_B)
    centered<-x-mean(x)+mu0
    out<-numeric(B)
    for(b in seq_len(B)) {
      z<-sample(centered,n,replace=TRUE); se<-sd(z)/sqrt(n)
      out[b]<-if(se>0) (mean(z)-mu0)/se else NA_real_
    }
    out
  },ignoreNULL=FALSE)
  ht_results<-reactive({
    x<-ht_sample();mu0<-input$ht_mu0;n<-length(x);se<-sd(x)/sqrt(n)
    observed<-if(se>0) (mean(x)-mu0)/se else NA_real_
    z<-ht_null(); z<-z[is.finite(z)]
    p<-if(length(z)>0 && is.finite(observed)) (1+sum(abs(z)>=abs(observed)))/(length(z)+1) else NA_real_
    list(t=observed,z=z,p=p,theory=if(is.finite(observed)) 2*pt(-abs(observed),df=n-1) else NA_real_)
  })
  output$ht_plot<-renderPlot({
    d<-ht_results();req(length(d$z)>10,is.finite(d$t))
    threshold<-abs(d$t)
    breaks<-sort(unique(c(pretty(range(c(d$z,-threshold,threshold)),n=35),-threshold,threshold)))
    h<-hist(d$z,breaks=breaks,plot=FALSE)
    cols<-rep(PAL$bootstrap,length(h$counts))
    cols[h$breaks[-1]<=-threshold | h$breaks[-length(h$breaks)]>=threshold]<-PAL$rose
    ymax<-max(h$density,dt(0,length(ht_sample())-1)) * 1.24
    plot(h,freq=FALSE,col=cols,border="black",las=1,bty="l",ylim=c(0,ymax),
         xlab="One-sample t-statistic",ylab="Density",main="Null-centred bootstrap test")
    curve(dt(x,df=length(ht_sample())-1),add=TRUE,col=PAL$blue,lwd=2)
    segments(c(-threshold,threshold),0,c(-threshold,threshold),ymax*.76,col=PAL$rose,lty=2,lwd=2)
    legend("topright",c("Bootstrap p-value tails","Observed test statistic","Theoretical t-distribution"),
           pch=c(15,NA,NA),col=c(PAL$rose,PAL$rose,PAL$blue),lty=c(0,2,1),lwd=c(0,2,2),
           pt.cex=1.1,bty="n",cex=.8)
  })
  output$ht_hypotheses <- renderUI({
    mu0 <- input$ht_mu0
    req(is.finite(mu0))
    div(class="result",
        strong("Hypotheses (two-sided test)"), tags$br(),
        mj(sprintf("\\[H_0: \\mu = %s \\qquad H_1: \\mu \\ne %s\\]",
                   fmt(mu0, 2), fmt(mu0, 2))))
  })
  output$ht_results<-renderUI({
    d<-ht_results();div(class="result",strong(mj(paste0("Observed \\(t\\)-statistic: ",fmt(d$t)))),tags$br(),
                        mj(paste0("Bootstrap \\(p\\)-value: ",fmt(d$p,4))),tags$br(),
                        mj(paste0("Theoretical \\(t\\)-test \\(p\\)-value: ",fmt(d$theory,4))),tags$br(),
                        sprintf("Number of finite null bootstrap statistics: %d",length(d$z)))
  })
  
  for(prefix in c("ex","co","cv")) local({pfx<-prefix
  observeEvent(input[[paste0(pfx,"_reset")]],{
    updateSelectInput(session,paste0(pfx,"_dist"),selected="lnorm")
    updateSelectInput(session,paste0(pfx,"_stat"),selected="mean")
    updateNumericInput(session,paste0(pfx,"_n"),value=30)
    updateSelectInput(session,paste0(pfx,"_B"),selected=2000)
    updateSliderInput(session,paste0(pfx,"_level"),value=.95)
    updateSliderInput(session,paste0(pfx,"_q"),value=.75)
    updateNumericInput(session,paste0(pfx,"_mu"),value=0)
    updateNumericInput(session,paste0(pfx,"_sigma"),value=1)
  })
  })
  observeEvent(input$ht_reset,{
    updateSelectInput(session,"ht_dist",selected="norm")
    updateNumericInput(session,"ht_n",value=30)
    updateSelectInput(session,"ht_B",selected=2000)
    updateNumericInput(session,"ht_mu0",value=0)
    updateNumericInput(session,"ht_mu",value=0)
    updateNumericInput(session,"ht_sigma",value=1)
  })
}

shinyApp(ui,server)
