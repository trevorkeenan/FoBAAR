REAL FUNCTION AssignLAI(lai,xfang,doy,lat,j,subDaily,par,PPFDsun,PPFDshd)
! Split canopy into sun and shade, and calculate intercepter radiation.

implicit none
	real, intent(in) :: lai,xfang,doy,lat,par
	real, intent(out) :: PPFDsun,PPFDshd
	integer, intent(in) :: j,subDaily
	
	integer :: nw
	
	real :: par2, albedo_par,parabs,parabs_laisun,parabs_laishade
	real :: parabs_per_laisun,parabs_per_laishade
	real :: cozen15,cozen45,cozen75,xk15,xk45,xk75,transd,extkd,radsol
	real :: sinbet,solext
	real :: tmprat,tmpr,tmpk,fdiff,fbeam,taul(2),rhol(2),rhos(2),pidiv,scatt(2),rhoch
	real :: rhoc(2,2),rhoc15,rhoc45,rhoc75,kpr(2,2),reff(2,2),Qabs(2,2),Qd0,Qb0,radabv(2)
	real :: sinlat,coslat,sindec,cosdec,A,B,coszen	
	real :: xphi1,xphi2,funG,extKb,pi
	
	par2 = 0
	albedo_par= 0
	parabs= 0
	parabs_laisun= 0
	parabs_laishade= 0
	parabs_per_laisun= 0
	parabs_per_laishade= 0
	cozen15= 0
	cozen45= 0
	cozen75= 0
	xk15= 0
	xk45= 0
	xk75= 0
	transd= 0
	extkd= 0
	radsol= 0
	sinbet= 0
	solext= 0
	tmprat= 0
	tmpr= 0
	tmpr= 0
	fdiff= 0
	fbeam= 0
	taul= 0
	rhol= 0
	rhos= 0
	pidiv= 0
	scatt= 0
	rhoch= 0
	rhoc(2,2)= 0
	rhoc15= 0
	rhoc45= 0
	rhoc75= 0
	kpr(2,2)= 0
	reff(2,2)= 0
	Qabs(2,2)= 0
	Qd0= 0
	Qb0= 0
	radabv= 0
	sinlat= 0
	coslat= 0
	sindec= 0
	cosdec= 0
	A= 0
	B= 0
	coszen= 0	
	xphi1= 0
	xphi2= 0
	funG= 0
	extKb= 0
	
	pi = 3.1416
	
	! 1 - calculate the solar zenith angle
	! calculations according to Goudriaan & van Laar 1994 p30
	      pidiv = pi/180.			! pi / 180
	
	! sine and cosine of latitude
	      sinlat = sin(lat)		! lat is in radians 
	      coslat = cos(lat)
	      
	! sine of maximum declination
	      sindec=-sin(23.45*pidiv)*cos(2.0*pi*(doy+10.0)/365.0)
	      cosdec=sqrt(1.-sindec*sindec)			! equation 3.5 Goudriaan & van Laar
	! terms A & B
	      A = sinlat*sindec
	      B = coslat*cosdec
	      !coszen = A+B*cos(2*pi*(j-(subDaily/2.))/(subDaily))   ! = solar zenith angle (equation 3.2 Goudriaan & van Laar)
          ! correction of Min Lee ('13)
          coszen = A+B*cos(2*pi*(((j-1)-(subDaily/2))/(subDaily))) ! = solar zenith angle (equation 3.2 Goudriaan & van Laar)

	
	
	! 2 - calculate the fraction of sunlit leaves
	! Ross-Goudriaan function for G(u) (see Wang et al. 2007 eq 13)
	xphi1 = 0.5 - 0.633*xfang -0.33*xfang*xfang
	xphi2 = 0.877 * (1.0 - 2.0*xphi1)
	funG=xphi1 + xphi2*coszen                   !G-function: Projection of unit leaf area in direction of beam
	
	if(coszen.gt.0) then                                 		! check if day or night
		extKb=funG/coszen                              		! beam extinction coeff - black leaves
		!AssignLAI = (1.0-exp(-extKb*lai))/extKb  ! LAI of sunlit leaves
        AssignLAI = (1.0-exp(-extKb*lai))/(extKb*lai)  ! LAI of sunlit leaves, Min Lee correction

		if (AssignLAI.gt.LAI)AssignLAI=LAI
		if (AssignLAI.lt.0)AssignLAI=0
	else
	  	extKb=100.
	  	AssignLAI=0.
	endif
	
	
	! Effective extinction coefficient for diffuse radiation Goudriaan & van Laar Eq 6.6)
	cozen15=cos(pidiv*15)
	cozen45=cos(pidiv*45)
	cozen75=cos(pidiv*75)
	xK15=xphi1/cozen15+xphi2
	xK45=xphi1/cozen45+xphi2
	xK75=xphi1/cozen75+xphi2
	transd=0.308*exp(-xK15*LAI)+0.514*exp(-xK45*LAI)+0.178*exp(-xK75*LAI)
	extkd=0.8 !(-1./LAI)*alog(transd)
	
	! Goudriaan theory as used in Leuning et al 1995 (Eq Nos from Goudriaan & van Laar, 1994)
	radsol=par
	sinbet=coszen
	solext=1370.0*(1.0+0.033*cos(2.0*pi*(doy-10.)/365.0))*sinbet
	tmprat=radsol/solext
	
	tmpR=0.847-1.61*sinbet+1.04*sinbet*sinbet
	tmpK=(1.47-tmpR)/1.66
	if(tmprat.le.0.22) fdiff=1.0
	if(tmprat.gt.0.22.and.tmprat.le.0.35) then
		fdiff=1.0-6.4*(tmprat-0.22)*(tmprat-0.22)
	endif
	
	if((tmprat.gt.0.35).and.(tmprat.le.tmpK)) then
		fdiff=1.47-1.66*tmprat
	endif
	
	if(tmprat.ge.tmpK) then
		fdiff=tmpR
	endif
	
	fbeam=1.0-fdiff
	tauL(1)=0.1                  ! leaf transmittance for vis
	tauL(2)=0.425                ! for NIR
	      
	rhoL(1)=0.1                  ! leaf reflectance for vis
	rhoL(2)=0.425                ! for NIR
	
 	rhoS(1)=0.1                  ! soil reflectance for vis
  	rhoS(2)=0.3                  ! for NIR - later function of soil water content
	
	do nw=1,2                                                         !nw:1=VIS, 2=NIR
		scatt(nw)=tauL(nw)+rhoL(nw)                                      !scattering coeff
		kpr(nw,1)=extKb*sqrt(1.-scatt(nw))                               !modified k beam scattered (6.20)
		kpr(nw,2)=extkd*sqrt(1.-scatt(nw))                               !modified k diffuse (6.20)
		rhoch=(1.-sqrt(1.-scatt(nw)))/(1.+sqrt(1.-scatt(nw)))            !canopy reflection black horizontal leaves (6.19)
		rhoc15=2.0*xK15*rhoch/(xK15+extkd)                                !canopy reflection (6.21) diffuse
		rhoc45=2.0*xK45*rhoch/(xK45+extkd)
		rhoc75=2.0*xK75*rhoch/(xK75+extkd)
		rhoc(nw,2)=0.057 !0.308*rhoc15+0.514*rhoc45+0.178*rhoc75, Min Lee correction
		rhoc(nw,1)=2.0*extKb/(extKb+extkd)*rhoch                          !canopy reflection (6.21) beam 
		reff(nw,1)=rhoc(nw,1)+(rhoS(nw)-rhoc(nw,1))*exp(-2.*kpr(nw,1)*LAI)                    !effective canopy-soil reflection coeff - beam (6.27)
		reff(nw,2)=rhoc(nw,2)+(rhoS(nw)-rhoc(nw,2)) *exp(-2.*kpr(nw,2)*LAI)                       !effective canopy-soil reflection coeff - diffuse (6.27)
	enddo
	      
	radabv(1) =par
	radabv(2) =0	! we only use par input ....
	do nw=1,1										!nw:1=VIS, 2=NIR
		Qd0=(1.-fbeam)*radabv(nw)                                          !diffuse incident radiation
		Qb0=fbeam*radabv(nw)                                               !beam incident radiation
!		Qabs(nw,2)=Qd0*(kpr(nw,2)*(1.-reff(nw,2))*exp(-kpr(nw,2)*LAI))+&   ! total absorbed radiation - shaded leaves, diffuse
!		&            Qb0*(kpr(nw,1)*(1.-reff(nw,1))*exp(-kpr(nw,1)*LAI)-&    !beam scattered
!		&             extKb*(1.-scatt(nw))*exp(-extKb*LAI))
		Qabs(nw,2)=Qd0*(kpr(nw,2)*(1.-reff(nw,2))*exp(-kpr(nw,2)*LAI))+&   ! total absorbed radiation - shaded leaves, diffuse
		&            Qb0*(kpr(nw,1)*(1.-reff(nw,1))*exp(-kpr(nw,1)*LAI)-&    !beam scattered
		&             extKb*(1.-scatt(nw))*exp(-extKb*LAI))
		Qabs(nw,1)=Qabs(nw,2)+extKb*Qb0*(1.-scatt(nw))                     ! total absorbed radiation - sunlit leaves 
	enddo
	
	PPFDsun=Qabs(1,1)
	PPFDshd=Qabs(1,2)
	
		!/* calculate PAR absorbed */
		! adopted from Biome-BGC version 4.2
		! method uses two variables, ext_coef and albedo
!		ext_coef = 0.7	
!		albedo = 0.1
!		k_par = ext_coef * 1.0;
!		albedo_par = albedo/3.0;
!		
!		Par2 = par * (1.0 - albedo_par);
!		parabs = par2 * (1.0 - exp(-k_par*lai));
!		
		!/* calculate the total PAR absorbed by the sunlit and
		!shaded canopy fractions */
!		parabs_laisun = k_par * par2 * Assignlaisun;
!		parabs_laishade = parabs - parabs_laisun;
!		
!		if (parabs_laishade.lt.0.0)then
!			parabs_laisun = parabs;
!			parabs_laishade = 0.0;
!		end if
!		
		! convert this to the PAR absorbed per unit LAI in the sunlit and 
		! shaded canopy fractions
!		if ((lai.gt.0.0).and.(Assignlaisun.gt.0.0))then
!			parabs_per_laisun = parabs_laisun/AssignLAISun;
!			parabs_per_laishade = parabs_laishade/(lai-AssignLAISun);
!		else
!			parabs_per_laisun = 0.
!			parabs_per_laishade = 0.
!		end if
!		
!		PPFDsun = parabs_per_laisun 
!		PPFDshd = parabs_per_laishade
!		
	! the following if radiation is read in in W/m2
	!/* calculate PPFD: assumes an average energy for PAR photon (EPAR, umol/J)
	!unit conversion: W/m2 --> umol/m2/s. */
	!ppfd_per_laisun = parabs_per_laisun *4.55 	! (umol/J) PAR photon energy ratio
	!ppfd_per_laishade = parabs_per_laishade *4.55

Return			
	
END

