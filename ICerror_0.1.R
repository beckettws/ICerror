as.ICdata = function(data, class = "iid"){
  if(class(data) == "paleoTS") class = "paleoTS"
  myICdata = list(data = data, class = class)

  class(myICdata) = c("ICdata",class)
  return(myICdata)
}



AIC.iid = function(myICfit, myICdata){
  k = length(myICfit$parameters) - length(myICfit$fix.arg)
  return(2*(k - logLik(myICfit, myICdata)))
}
AIC.ICpaleo = function(myICfit, myICdata){
  k = myICfit$data$K
  return(2*(k - logLik(myICfit, myICdata)))
}

AICc = function(myICfit, myICdata){
  UseMethod("AICc", myICfit)
}
AICc.ICpaleo = function(myICfit, myICdata){
  myICfit$data$AICc
}
logLik.ICfit = function(myICfit, myICdata) NextMethod("logLik", myICfit)
logLik.iid = function(myICfit, myICdata){
  NextMethod("logLik", myICfit)
}
logLik.norm  =  function(myICfit, myICdata){
  R = dnorm(myICdata$data, myICfit$parameters["mean"], myICfit$parameters["sd"])
  return(sum(log(R)))
}
logLik.exp  =  function(myICdata, myICfit){
  R = dnorm(myICdata$data, myICfit$parameters["mean"], myICfit$parameters["sd"])
  return(sum(log(R)))
}
logLik.ICpaleo = function(myICfit, myICdata){
  myICfit$data$logL
}

as.ICdist = function(name = "norm", fix.arg = list()){

  iidlist = c("norm", "exp")
  paleolist = c("GRW", "URW",
                   "Stasis", "StrictStasis",
                   "Punc", "StasisURW", "StasisGRW",
                   "URWStasis", "GRWStasis")

  class = ifelse(name %in% paleolist, "ICpaleo", "iid")
  class = gsub("-", "", class)
  fix = ""
  if (length(fix.arg) > 0) {
    fix = paste(": ",names(fix.arg), " fixed", sep = "")}
  desc = paste(name, fix, sep = "")
  mydist = list(name = name, fix.arg = fix.arg, desc = desc)

  class(mydist) = c("ICdist", class, name)
  return(mydist)
}

as.ICfit = function(name = "norm", parameters = c("mean" = 0, "sd" = 1), fix.arg = list(), data = list()){
  iidlist = c("norm", "exp")
  paleolist = c("GRW", "URW",
                "Stasis", "StrictStasis",
                "Punc", "Punc-1" ,"Stasis-URW", "Stasis-GRW",
                "URW-Stasis", "GRW-Stasis")

 name =  ifelse(name == "Punc-1", "Punc", name)
  class = ifelse(name %in% paleolist, "ICpaleo", "iid")
  ifelse(class == "Punc-1", "Punc", class)
  name = sub("-", "", name)
  myfit = list(name = name, parameters = parameters, fix.arg = fix.arg, data = data)
  class(myfit) = c("ICfit", class, name)
  return(myfit)
}

ICpaleofit = function(paleoTSfit){
  as.ICfit(paleoTSfit$modelName, paleoTSfit$parameters, data = paleoTSfit)
}
dists =  function(distlist) unlist(lapply(distlist, function(x) tail(class(x),1)))

fit = function(ICdata, ICdist, ...){
  UseMethod("fit", ICdist)
}

  fit.ICdist = function(ICdata, ICdist, ...) NextMethod("fit", ICdist)
    fit.iid = function(ICdata, ICdist, ...){
      fit = fitdist(ICdata$data, ICdist$name, fix.arg = ICdist$fix.arg)
      param = c(list(fit$estimate), fit$fix.arg, recursive = TRUE)
      myICfit = as.ICfit(ICdist$name, param, fit$fix.arg)
      # class(myICfit) = append(class(myICfit), ICdist$name)
      return(myICfit)
    }



  fit.ICpaleo = function(ICdata, ICdist, ...){
   NextMethod("fit", ICdist)
     # fitlist = list(
    #   "GRW" = fit.GRW,
    #   "URW" = opt.URW,
    #   "Stasis" = fit.Stasis
    #   "StrictStasis" = fit.StrictStasis
    # )
  }

  fit.GRW = function(ICdata, ICdist,...){
    ICpaleofit(opt.GRW(ICdata$data, pool = FALSE))
  }
  fit.URW = function(ICdata, ICdist,...){
    ICpaleofit( opt.URW(ICdata$data, pool = FALSE))
  }
  fit.Stasis = function(ICdata, ICdist,...){
    ICpaleofit(opt.Stasis(ICdata$data, pool = FALSE))
  }
  fit.StrictStasis = function(ICdata, ICdist,...){
    ICpaleofit(opt.StrictStasis(ICdata$data, pool = FALSE))
  }
  fit.Punc = function(ICdata, ICdist,...){
    ICpaleofit(fitGpunc(ICdata$data, ng=2, pool = FALSE, silent = TRUE))
  }
  fit.StasisURW = function(ICdata, ICdist,...){
    ICpaleofit(fitModeShift(ICdata$data, order="Stasis-RW", rw.model="URW", silent = TRUE, pool = FALSE))
  }
  fit.StasisGRW = function(ICdata, ICdist,...){
    ICpaleofit( fitModeShift(ICdata$data, order="Stasis-RW", rw.model="GRW",silent = TRUE, pool = FALSE))
  }
  fit.URWStasis = function(ICdata, ICdist,...){
    ICpaleofit(fitModeShift(ICdata$data, order="RW-Stasis", rw.model="URW", silent = TRUE, pool = FALSE))
  }
  fit.GRWStasis = function(ICdata, ICdist,...){
    ICpaleofit(fitModeShift(ICdata$data, order="RW-Stasis", rw.model="GRW", silent = TRUE,pool = FALSE))
  }

score = function(ICfit, ICdata, IC = "AIC"){
  f = match.fun(IC)
  f(ICfit, ICdata)
}

#  score.paleoTS = function(ICdata, ICdist, IC = "AIC"){
    UseMethod(IC, iCdist)
  }

scores = function(ICdata, ICdistlist, IC = "AIC"){
  lapply(iCdistlist, function(x) score(data, x, IC))
}

NPBoot = function(ICdata, N = 1000){
  UseMethod("NPBoot", ICdata)
}
  NPBoot.ICdata = function(ICdata, N = 1000){
    NextMethod("NPBoot", ICdata)
}
  NPBoot.iid = function(ICdata, N = 1000){

      lapply(1:N, function(x) as.ICdata(sample(ICdata$data, replace = TRUE), class = "iid"))
}


  NPBoot.paleoTS = function(my.ICpaleo, N = 1000){
    samples = rnorm(sum(my.ICpaleo$data$nn) * N, rep(rep(my.ICpaleo$data$mm, my.ICpaleo$data$nn), N), rep(rep(my.ICpaleo$data$vv, my.ICpaleo$data$nn), N))
    samples.split1 = split(samples, ceiling(seq_along(samples)/sum(my.ICpaleo$data$nn)))
    samples.split2 = lapply(samples.split1, function(x) splitAt(x, pos = 1 + cumsum(my.ICpaleo$data$nn)))
    mmvv = lapply(samples.split2, function(x) cbind("mm" = as.numeric(lapply(x, mean)), "vv" = as.numeric(lapply(x, var))))
    return(lapply(mmvv, function(x) as.ICdata(as.paleoTS(x[,1], x[,2], my.ICpaleo$data$nn, my.ICpaleo$data$tt ))))
  }

newFits = function(ICdatalist, ICdistlist){
  UseMethod(ICdistlist[[1]])
}
  newFits.iid = function(ICdatalist, ICdistlist){

  }
  newFits.paleoTS = function(ICdatalist, ICdistlist){

  }


  splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))


samp = function(n, ICdist){
  UseMethod("samp", ICdist)
}

  samp.iid = function(n, ICdist){
    NextMethod("samp", ICdist)
    }
    samp.norm = function(n, ICdist){
      rnorm(n, ICdist$parameters["mean"], ICdist$parameters["sd"])
    }
    samp.exp = function(n, ICdist){
      rexp(n, ICdist$parameters["rate"])
    }

  samp.ICpaleo = function(n, ICdist){
    NextMethod("samp", ICdist)
  }
    samp.GRW = function(n, ICdist){
      sim.GRW(ns = 20, ms = ICdist$parameters["mstep"], vs = ICdist$parameters["vstep"])
    }
    samp.URW = function(n, ICdist){
      sim.GRW(ns = 20, ms = 0, vs = ICdist$parameters["vstep"])
    }
    samp.Stasis = function(n, ICdist){
      sim.Stasis(20,ICdist$parameters["theta"],ICdist$parameters["omega"], 1)
    }
    samp.StrictStasis = function(n, ICdist){
      sim.Stasis(20,ICdist$parameters["theta"],0, 1)
    }
    samp.Punc = function(n, ICdist){
      theta = ICdist$parameters[grepl("theta",names(ICdist$parameters))]
      sim.punc(ns = c(10, 10),theta, omega = rep(ICdist$parameters["omega"],length(theta)),
               nn = rep(30, sum(ns)), tt = 0:(sum(ns)-1), vp = 1)
    }


    sim.punc<- function (ns=c(10,10), theta=c(0,1), omega=rep(0,length(theta)), nn=rep(30,sum(ns)), tt=0:(sum(ns)-1), vp=1)
      # simulate punctuated sequence; theta and omega are vectors of paramters
      # ns is vector of ns in each sub-sequence
    {
      nr<- length(theta)
      xl<- list()
      for (i in 1:nr)
      {
        if (i==1)
        {  start.i<- 1
        end.i<- ns[i]  }
        else
        {
          start.i<- sum(ns[1:(i-1)])+1
          end.i<- start.i + ns[i] -1
        }

        xl[[i]]<- sim.Stasis(ns=ns[i], theta=theta[i], omega=omega[i], vp=vp, nn=nn[start.i:end.i], tt[start.i:end.i])
      }

      y<- cat.paleoTS(xl)
      y$label<- "Created by sim.punc()"
      shft<- cumsum(ns[-nr])+1
      y$genpars<- c(theta, omega, shft)
      names(y$genpars)<- c(paste("theta",1:nr,sep=""), paste("omega",1:nr,sep=""), paste("shift",1:(nr-1),sep=""))
      return (y)
    }

    ### How to do these???
    samp.StasisURW = function(n, ICdist){
      ns1 = ICdist$parameters["shift1"]-1
      ns2 = ICdist$data$n - ns1
      ns = c(ns1, ns2,0)
      tt2 = (ns1+1):sum(ns1,ns2)
      theta = ICdist$parameters["theta"]
      omega = ICdist$parameters["omega"]
      ms =  0
      vs = ICdist$parameters["vstep"]
      p1 = sim.Stasis(ns1, theta, omega)
      p2 = sim.GRW(ns2, ms, vs, tt = tt2)
      p2$MM =  p2$MM + p1$MM[ns1]
      cat.paleoTS(list(p1, p2))
    }

    samp.StasisGRW = function(n,ICdist){
      ns1 = ICdist$parameters["shift1"]-1
      ns2 = ICdist$data$n - ns1
      ns = c(ns1, ns2,0)
      tt2 = (ns1+1):sum(ns1,ns2)
      theta = ICdist$parameters["theta"]
      omega = ICdist$parameters["omega"]
      ms =  ICdist$parameters["mstep"]
      vs = ICdist$parameters["vstep"]
      p1 = sim.Stasis(ns1, theta, omega)
      p2 = sim.GRW(ns2, ms, vs, tt = tt2)
      p2$MM =  p2$MM + p1$MM[ns1]
      cat.paleoTS(list(p1, p2))
    }

    samp.URWStasis = function(n, ICdist){
      ns1 = ICdist$parameters["shift1"]-1
      ns2 = ICdist$data$n - ns1
      ns = c(ns1, ns2)
      tt2 = (ns1+1):sum(ns1,ns2)
      theta = ICdist$parameters["theta"]
      omega = ICdist$parameters["omega"]
      ms =  0
      vs = ICdist$parameters["vstep"]
      p1 = sim.GRW(ns1, ms, vs)
      p2 = sim.Stasis(ns2, theta, omega, tt = tt2)
      p2$MM =  p2$MM + p1$MM[ns1]
      cat.paleoTS(list(p1, p2))

    }
    samp.GRWStasis = function(n, ICdist){
      ns1 = ICdist$parameters["shift1"]-1
      ns2 = ICdist$data$n - ns1
      ns = c(ns1, ns2)
      tt2 = (ns1+1):sum(ns1,ns2)
      theta = ICdist$parameters["theta"]
      omega = ICdist$parameters["omega"]
      ms =  ICdist$parameters["mstep"]
      vs = ICdist$parameters["vstep"]
      p1 = sim.GRW(ns1, ms, vs)
      p2 = sim.Stasis(ns2, theta, omega, tt = tt2)
      p2$MM =  p2$MM + p1$MM[ns1]
      cat.paleoTS(list(p1, p2))

    }


    sim.sg = function(ns=c(20,20), theta=0, omega=1, ms=1, vs=0.1, nn=rep(30, sum(ns)), tt=0:(sum(ns)-1), vp=1)


    sim.sgs2 <- function (ns=c(20,20,20), theta=0, omega=1, ms=1, vs=0.1, nn=rep(30, sum(ns)), tt=0:(sum(ns)-1), vp=1)
      # simulate stasis-grw-stasis sequence, take theta2 to be final value after grw part
    {
      xl<- list()
      for (i in 1:3)
      {
        if (i==1)
        {  start.i<- 1
        end.i<- ns[i]
        }
        else
        {
          start.i<- sum(ns[1:(i-1)])+1
          end.i<- start.i + ns[i] -1
        }

        if (i==2)
          xl[[i]] = 0
          if (start.i >= end.i)
            xl[[i]] =  sim.GRW(ns[2], ms, vs, nn=nn[start.i:end.i], tt=tt[start.i:end.i], vp=vp)


        else
          xl[[i]] = 0
        if (start.i >= end.i)
          xl[[i]]<- sim.Stasis(ns=ns[i], theta=theta, omega=omega, vp=vp, nn=nn[start.i:end.i], tt[start.i:end.i])
      }

      ## add offsets
      xl[[2]]$mm<- xl[[2]]$mm + xl[[1]]$MM[ns[1]]
      xl[[3]]$mm<- xl[[3]]$mm + xl[[2]]$MM[ns[2]]

      y<- cat.paleoTS(xl)
      y$label<- "Created by sim.sgs()"
      y$genpars <- c(theta, omega, ms, vs)
      names(y$genpars)<- c("theta","omega", "ms","vs")
      return(y)
    }

ICscores = function(ICdata, ICdistlist, IC = "AIC", N = 10){
  ###Fit dists to data

  datatype = class(ICdata)[2]
  distnames = sapply(ICdistlist, function(x) x$desc)
  obs.fits = lapply(ICdistlist, function(x) fit(ICdata, x))

  ###Observed Scores
  obs.scores = sapply(obs.fits, function(x) score(x, ICdata))
  names(obs.scores) = distnames

  ###NP Bootstrap

  bootdata = NPBoot(ICdata, N)

  ###Fit New Models

  bootfits = lapply(ICdistlist, function(x) lapply(bootdata, function(y )fit(y, x)))
  names(bootfits) = distnames

  ###Parametric Resample
  resamples = lapply(bootfits, function(x) lapply(x, function(y) as.ICdata(samp(length(ICdata$data), y), class = datatype)))
  names(resamples) = distnames
  ###New Scores
  bootscores = lapply(resamples, function(x) t(sapply(x, function(y) sapply(ICdistlist, function(z) score(fit(y, z), y)))))
  bootscores = lapply(bootscores, function(x) {
    colnames(x) = distnames
     x})
  names(bootscores) = distnames


  names(bootscores) = paste0(names(bootscores), rep(" true"))
  return(bootscores)

}
