# xy_str <- function(e) {
#   if(is.null(e)) return("NULL\n")
#   paste0("x=", round(e$x, 1), " y=", round(e$y, 1), "\n")
# }

###### response time and metrics

m_t_complete <- function( t , k1 , k2 , k3 ) {
  if( k2 != k3 ) {
    k1 / k3 * ( k2 * ( 2 - exp( - k3 * t ) ) -  k3 * ( 2 - exp( - k2 * t ) ) ) / ( k2 - k3 )
  } else {
    k1 / k3 * ( 2 - exp( - k3 * t )*( 1 + k3 * t ) )
  }
}

m_t_simple <- function( t , k1 , k3 ) {
  k1 / k3 * ( 2 - exp( - k3 * t ) )
}

hrtime <- function( k2 , k3 , maxerr=1e-2) 
{
  errfunHRfullsys <- function(t,k2,k3) 
    if( k2 != k3 ) {
      ( -k2 * exp(-k3*t) + k3 * exp(-k2 * t) + 1/2*(k2-k3) ) ^ 2
    } else {
      ( 1 - 2*exp( -k2 * t )*( k2 * t + 1 ) )^2
    }

  HRnoprocessing <- log(2) / k3 
  HRfullsys <- optimize(errfunHRfullsys, c(0,100*HRnoprocessing), 
    k2 = k2 , k3 = k3 )

  relerror <- HRfullsys$objective / k2 
  if( relerror < maxerr )
    return( HRfullsys$minimum )
  else NA
  
}

tau_fun <- function( k2 , k3 , maxerr=1e-2) {
  if( is.finite(k2) & is.finite(k3) ) {
    HRnoprocessing <- log(2) / k3 
    HRfullsys <- hrtime( k2 , k3 , maxerr )
    return( HRfullsys / HRnoprocessing  )
  } else {
    return( NaN )
  }
}

delta_fun <- function( k1 , k2 , k3 )
  if( k2 != k3 ) {
    k1 * ( 1 - 2^( 1 - k2 / k3 ) ) / ( 2 * ( k2 - k3 ) )
  } else {
    k1 *log(2) / 2 / k2  
  }


########################################################
## plot the dynamics of mRNA following an induction ####
## of the synthesis rate of 2 fold, with or without ####
## processing ##########################################
########################################################

shinyProcessingDelayPlot <- function(k1, k2, k3, absval, metrics) {

  simulation_time <- seq(0,1,length.out=1000)

  ## make simulations with and without processing
  
  mature_complete <- sapply(simulation_time, function(t) 
    m_t_complete(t, k1, k2, k3) )
  mature_simple <- sapply(simulation_time, function(t) 
    m_t_simple(t, k1, k3) )

  # normalize profiles relative to the inital state

  mature_complete_relative <- mature_complete/mature_complete[1]
  mature_simple_relative <- mature_simple/mature_simple[1]

  # calculate the half response time

  half_resp_time <- hrtime(k2, k3)
  half_resp_time_noproc <- log(2) / k3

  # calculate mRNA levels at half response times
  
  mat_levels <- m_t_complete(half_resp_time_noproc, k1, k2, k3)
  mat_levels_ht <- m_t_complete(half_resp_time, k1, k2, k3)
  mat_levels_noproc <- m_t_simple(half_resp_time_noproc, k1, k3)
  
  rel_mat_levels <- mat_levels/mature_complete[1]
  rel_mat_levels_ht <- mat_levels_ht/mature_complete[1]
  rel_mat_levels_noproc <- mat_levels_noproc/mature_simple[1]
  
  # calculate metric 1 and metric 2
  
  # tau_value <- half_resp_time/half_resp_time_noproc
  tau_value <- tau_fun(k2, k3)
  # delta_value <- mat_levels_noproc - mat_levels
  delta_value <- delta_fun(k1, k2, k3)
  
  if( absval ) {
    labeltag <- ''
  } else {
    labeltag <- 'relative'
    mature_complete <- mature_complete_relative
    mature_simple <- mature_simple_relative
    mat_levels <- rel_mat_levels
    mat_levels_ht <- rel_mat_levels_ht
    mat_levels_noproc <- rel_mat_levels_noproc
  }
  
  # start plot routine

  par(mar=c(5,4,4,2)+.1)
  matplot(simulation_time, 
          cbind(mature_complete, mature_simple),
          ylim=c(min(mature_simple), max(mature_simple)*1.05),
          xaxs='i', yaxs='i', xaxt='n',
          ylab = paste(labeltag, 
            'mature RNA upon doubling of synthesis'), 
          type='l', xlab='', col=1, lwd=2)
    
  if( metrics ) {
    
    # draw metric 1 and metric 2 segments on the plot
    
    segments(half_resp_time_noproc,mat_levels,half_resp_time_noproc,mat_levels_noproc, 
             col='darkgoldenrod1', lty=1, lwd=3)
    
    # draw metric 1 and metric 2 labels on the plot
    
    text(half_resp_time_noproc+.028,(mat_levels+mat_levels_noproc)/2, 
         expression(Delta), col='darkgoldenrod1', cex=1.5)

    title( bquote(tau ~ '=' ~ frac("t"[1/2] ~ "'", "t"[1/2]) ~ '=' ~ 
      .(round(tau_value,2)) ~ ';' ~ Delta ~ '=' ~ .(round(delta_value,2)) ) )

    # draw dotted lines at the two half response times (with and without processing)
    
    segments(0,mat_levels_noproc, half_resp_time,mat_levels_noproc,
            col='black', lty=3, lwd=1)  
    segments(half_resp_time_noproc,0,half_resp_time_noproc,mat_levels_noproc,
            col='black', lty=3, lwd=1)
    segments(half_resp_time,0,half_resp_time,mat_levels_ht,
            col='black', lty=3, lwd=1)

    # draw t1/2 signs on the x-axis
    
    at <- c(half_resp_time_noproc,half_resp_time, seq(0,1,by=.2))
    labels <- c(expression("t"[1/2]),expression("t"[1/2] ~ "'"), as.character(seq(0,.8,by=.2)), '1 hour')
    axis(1,at=at,labels=labels,las=2)

    
  } else {

    title( bquote(tau ~ '=' ~ .(round(tau_value,2)) ~ ';' ~ 
      Delta ~ '=' ~ .(round(delta_value,2)) ) )

    # draw t1/2 signs on the x-axis
  
    at <- seq(0,1,by=.2)
    labels <- c(as.character(seq(0,.8,by=.2)), '1 hour')
    axis(1,at=at,labels=labels,las=2)

  }
  
  # legend

  legend('bottomright', lty=2:1, bty='n', legend=c('istantaneous',
    paste0('measured (', round(k2,2),')')), title='processing:')
  
  return(list(tau_value=tau_value, delta_value=delta_value))

}

