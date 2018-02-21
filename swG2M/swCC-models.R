swcc.par <- list(
  hNCC = list(ic = c('z0'=1.8, 'z2'=1.2, 'y0'=1.2, 'y2'=1.8, 'q0'=1.2, 'q2'=1.8, 'p0'=1.8, 'p2'=1.2),
                 pp = c("nt"=1, "nd"=1, "hph"=0.5, "hkn"=0.5,
                        "Z"=3.0, "Y"=3.0, "Q"=3.0, "P"=3.0,
                        "k1"=1.0)),
  NCC = list(ic = c('z0'=1.8, 'z2'=1.2, 'y0'=1.2, 'y2'=1.8, 'q0'=1.2, 'q2'=1.8, 'p0'=1.8, 'p2'=1.2, 'r0'=1.8, 'r2'=1.2, 's0'=1.2, 's2'=1.8),
                 pp = c("nt"=1, "nd"=1, "hph"=0.5, "hkn"=0.5,
                        "Z"=3.0, "Y"=3.0, "Q"=3.0, "P"=3.0, "R"=3.0, "S"=3.0,
                        "k1"=1.0)),
  NCCchk = list(ic = c('z0'=1.8, 'z2'=1.2, 'y0'=1.2, 'y2'=1.8, 'q0'=1.2, 'q2'=1.8, 'p0'=1.8, 'p2'=1.2, 'r0'=1.8, 'r2'=1.2, 's0'=1.2, 's2'=1.8),
                  pp = c("nt"=8, "nd"=1, "ch"=0.5, "hph"=0.05,"hkn"=0.05,
                         "Z"=3.0, "Y"=3.0, "Q"=3.0, "P"=3.0, "R"=3.0, "S"=3.0,
                         "k1"=1.0))
)


swcc.mod <- list(
  
  hNCC = function(tt,x,k){
    with(as.list(c(x,k)),{
      dz0 = -k1*nd*y0*z0 + k1*nt*z0*(Z-z0-z2) - k1*hph*z0 + k1*hkn*(Z-z0-z2)
      dz2 = k1*nd*y0*(Z-z0-z2) - k1*nt*z0*z2 + k1*hph*(Z-z0-z2) - k1*hkn*z2
      dq0 = k1*nd*y0*(Q-q0-q2) - k1*nt*z0*q0 + k1*hph*(Q-q0-q2) - k1*hkn*q0
      dq2 = -k1*nd*y0*q2 + k1*nt*z0*(Q-q0-q2) - k1*hph*q2 + k1*hkn*(Q-q0-q2)
      dp0 = k1*nt*z0*(P-p0-p2) - k1*nd*y0*p0 - k1*hph*p0 + k1*hkn*(P-p0-p2)
      dp2 = -k1*nt*z0*p2 + k1*nd*y0*(P-p0-p2) + k1*hph*(P-p0-p2) - k1*hkn*p2
      dy0 = k1*nd*q0*(Y-y0-y2) - k1*nt*p0*y0 + k1*hph*(Y-y0-y2) - k1*hkn*y0
      dy2 = -k1*nd*q0*y2 + k1*nt*p0*(Y-y0-y2) - k1*hph*y2 + k1*hkn*(Y-y0-y2)
      return(list(c(dz0,dz2,dy0,dy2,dq0,dq2,dp0,dp2)))
    })},
  
  rNCC = function(tt,x,k){
    with(as.list(c(x,k)),{
      dz0 = -k1*nt*s0*z0 + k1*nd*r0*(Z-z0-z2) + k1*hkn*(Z-z0-z2) - k1*hph*z0
      dz2 = k1*nt*s0*(Z-z0-z2) - k1*nd*r0*z2 - k1*hkn*z2 + k1*hph*(Z-z0-z2)
      ds0 = k1*nd*y0*(S-s0-s2) - k1*nt*z0*s0 + k1*hkn*(S-s0-s2) - k1*hph*s0 
      ds2 = -k1*nd*y0*s2 + k1*nt*z0*(S-s0-s2) - k1*hkn*s2 + k1*hph*(S-s0-s2) 
      dr0 = k1*nt*z0*(R-r0-r2) - k1*nd*y0*r0 + k1*hph*(R-r0-r2) - k1*hkn*r0
      dr2 = -k1*nt*z0*r2 + k1*nd*y0*(R-r0-r2) - k1*hph*r2 + k1*hkn*(R-r0-r2)
      dq0 = k1*nd*y0*(Q-q0-q2) - k1*nt*z0*q0 - k1*hkn*q0 + k1*hph*(Q-q0-q2)
      dq2 = -k1*nd*y0*q2 + k1*nt*z0*(Q-q0-q2) + k1*hkn*(Q-q0-q2) - k1*hph*q2
      dp0 = k1*nt*z0*(P-p0-p2) - k1*nd*y0*p0 + k1*hkn*(P-p0-p2) - k1*hph*p0
      dp2 = -k1*nt*z0*p2 + k1*nd*y0*(P-p0-p2) - k1*hkn*p2 + k1*hph*(P-p0-p2)
      dy0 = k1*nd*q0*(Y-y0-y2) - k1*nt*p0*y0 - k1*hkn*y0 + k1*hph*(Y-y0-y2)
      dy2 = -k1*nd*q0*y2 + k1*nt*p0*(Y-y0-y2) + k1*hkn*(Y-y0-y2) - k1*hph*y2
      return(list(c(dz0,dz2,dy0,dy2,dq0,dq2,dp0,dp2,dr0,dr2,ds0,ds2)))
    })},
  
  aNCC = function(tt,x,k){
    with(as.list(c(x,k)),{
      dz0 = -k1*nd*s0*z0 + k1*nt*r0*(Z-z0-z2) + k1*hkn*(Z-z0-z2) - k1*hph*z0
      dz2 = k1*nd*s0*(Z-z0-z2) - k1*nt*r0*z2 - k1*hkn*z2 + k1*hph*(Z-z0-z2)
      ds0 = k1*nd*y0*(S-s0-s2) - k1*nt*z0*s0 + k1*hkn*(S-s0-s2) - k1*hph*s0 
      ds2 = -k1*nd*y0*s2 + k1*nt*z0*(S-s0-s2) - k1*hkn*s2 + k1*hph*(S-s0-s2) 
      dr0 = k1*nt*z0*(R-r0-r2) - k1*nd*y0*r0 + k1*hph*(R-r0-r2) - k1*hkn*r0
      dr2 = -k1*nt*z0*r2 + k1*nd*y0*(R-r0-r2) - k1*hph*r2 + k1*hkn*(R-r0-r2)
      dq0 = k1*nd*y0*(Q-q0-q2) - k1*nt*z0*q0 - k1*hkn*q0 + k1*hph*(Q-q0-q2)
      dq2 = -k1*nd*y0*q2 + k1*nt*z0*(Q-q0-q2) + k1*hkn*(Q-q0-q2) - k1*hph*q2
      dp0 = k1*nt*z0*(P-p0-p2) - k1*nd*y0*p0 + k1*hkn*(P-p0-p2) - k1*hph*p0
      dp2 = -k1*nt*z0*p2 + k1*nd*y0*(P-p0-p2) - k1*hkn*p2 + k1*hph*(P-p0-p2)
      dy0 = k1*nd*q0*(Y-y0-y2) - k1*nt*p0*y0 - k1*hkn*y0 + k1*hph*(Y-y0-y2)
      dy2 = -k1*nd*q0*y2 + k1*nt*p0*(Y-y0-y2) + k1*hkn*(Y-y0-y2) - k1*hph*y2
      return(list(c(dz0,dz2,dy0,dy2,dq0,dq2,dp0,dp2,dr0,dr2,ds0,ds2)))
    })},
  
    chkg2m = function(tt,x,k){
    with(as.list(c(x,k)),{
      dz0 = -k1*nt*s0*z0 + k1*nd*r0*(Z-z0-z2) + k1*hkn*(Z-z0-z2) - k1*hph*z0
      dz2 = k1*nt*s0*(Z-z0-z2) - k1*nd*r0*z2 - k1*hkn*z2 + k1*hph*(Z-z0-z2)
      ds0 = k1*nd*y0*(S-s0-s2) - k1*nt*z0*s0 + k1*hkn*(S-s0-s2) - k1*hph*s0 + k1*nt*ch*(S-s0-s2)
      ds2 = -k1*nd*y0*s2 + k1*nt*z0*(S-s0-s2) - k1*hkn*s2 + k1*hph*(S-s0-s2) - k1*nt*ch*s2
      dr0 = k1*nt*z0*(R-r0-r2) - k1*nd*y0*r0 + k1*hph*(R-r0-r2) - k1*hkn*r0 - k1*nt*ch*r0
      dr2 = -k1*nt*z0*r2 + k1*nd*y0*(R-r0-r2) - k1*hph*r2 + k1*hkn*(R-r0-r2) + k1*nt*ch*(R-r0-r2)
      dq0 = k1*nd*y0*(Q-q0-q2) - k1*nt*z0*q0 - k1*hkn*q0 + k1*hph*(Q-q0-q2)
      dq2 = -k1*nd*y0*q2 + k1*nt*z0*(Q-q0-q2) + k1*hkn*(Q-q0-q2) - k1*hph*q2
      dp0 = k1*nt*z0*(P-p0-p2) - k1*nd*y0*p0 + k1*hkn*(P-p0-p2) - k1*hph*p0
      dp2 = -k1*nt*z0*p2 + k1*nd*y0*(P-p0-p2) - k1*hkn*p2 + k1*hph*(P-p0-p2)
      dy0 = k1*nd*q0*(Y-y0-y2) - k1*nt*p0*y0 - k1*hkn*y0 + k1*hph*(Y-y0-y2)
      dy2 = -k1*nd*q0*y2 + k1*nt*p0*(Y-y0-y2) + k1*hkn*(Y-y0-y2) - k1*hph*y2
      return(list(c(dz0,dz2,dy0,dy2,dq0,dq2,dp0,dp2,dr0,dr2,ds0,ds2)))
    })},
  phogkp = function(tt,x,k){
    with(as.list(c(x,k)),{
      dz0 = -k1*nd*s0*z0 + k1*nt*r0*(Z-z0-z2) + k1*hkn*(Z-z0-z2) - k1*hph*z0
      dz2 = k1*nd*s0*(Z-z0-z2) - k1*nt*r0*z2 - k1*hkn*z2 + k1*hph*(Z-z0-z2)
      ds0 = k1*nd*y0*(S-s0-s2) - k1*nt*z0*s0 + k1*hkn*(S-s0-s2) - k1*hph*s0 + k1*nd*ch*(S-s0-s2)
      ds2 = -k1*nd*y0*s2 + k1*nt*z0*(S-s0-s2) - k1*hkn*s2 + k1*hph*(S-s0-s2) - k1*nd*ch*s2
      dr0 = k1*nt*z0*(R-r0-r2) - k1*nd*y0*r0 + k1*hph*(R-r0-r2) - k1*hkn*r0 - k1*nd*ch*r0
      dr2 = -k1*nt*z0*r2 + k1*nd*y0*(R-r0-r2) - k1*hph*r2 + k1*hkn*(R-r0-r2) + k1*nd*ch*(R-r0-r2)
      dq0 = k1*nd*y0*(Q-q0-q2) - k1*nt*z0*q0 - k1*hkn*q0 + k1*hph*(Q-q0-q2)
      dq2 = -k1*nd*y0*q2 + k1*nt*z0*(Q-q0-q2) + k1*hkn*(Q-q0-q2) - k1*hph*q2
      dp0 = k1*nt*z0*(P-p0-p2) - k1*nd*y0*p0 + k1*hkn*(P-p0-p2) - k1*hph*p0
      dp2 = -k1*nt*z0*p2 + k1*nd*y0*(P-p0-p2) - k1*hkn*p2 + k1*hph*(P-p0-p2)
      dy0 = k1*nd*q0*(Y-y0-y2) - k1*nt*p0*y0 - k1*hkn*y0 + k1*hph*(Y-y0-y2)
      dy2 = -k1*nd*q0*y2 + k1*nt*p0*(Y-y0-y2) + k1*hkn*(Y-y0-y2) - k1*hph*y2
      return(list(c(dz0,dz2,dy0,dy2,dq0,dq2,dp0,dp2,dr0,dr2,ds0,ds2)))
    })}
  )
