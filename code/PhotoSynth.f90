REAL FUNCTION PhotoSynth(airT,ppfdACT,VPD,Ca,Nit,totLAI,Vcmax,EaVcmax,EdVcmax,EaJmax,EdJmax,SJmax,Rd,Rdt,VQ10,gsD0,tetaph,g1,Gc) !,tetaph,gsD0

implicit none
  	real, intent(in) :: airT,ppfdAct,VPD,Ca,Nit,totLAI
	real, intent(in) :: Vcmax,EaVcmax,EdVcmax,EaJmax,EdJmax,SJmax,Rd,VQ10    	! Parameters to be optomised
	real, intent(in) :: g1,gsD0,tetaph
	real, intent(out) :: Rdt,Gc
	real :: Gammast25=42.22		! values taken from Bernacchi et al. 2001
	real :: EaGammast=37830		! ""
	real :: Kc=404.9				! ""
	real :: Ko=278400				! ""
	real :: EaKo=36380		! ""
	real :: EaKc=79430		! ""
	!real :: tetaph=0.5		! curvature AN/ppfd
	real :: pabs = 0.85;    	!/* (DIM) fPAR effectively absorbed by PSII */
	real :: ppe=2.6			!/* (mol/mol) photons absorbed by PSII per e- transported */
     	!kc 58,520 (2.21), ko 58,520 (2.21), Kc 59,356 (2.24), Ko 35,948 (1.63)
		!	real :: Vcmax=75
		!	real :: EaVcmax = 51560
		!	real :: EdVcmax=226000
		!	real :: EaJmax=43790
		!	real :: EdJmax=200000
		!	real :: SJmax=710
		!	real :: Rd=0.75
		!	real :: VQ10=1.2
	 
	real :: oi=207305		! umoles/mol (ppm)
	real :: gs0=0.01
	!real :: g1 = 11
	!real :: gsD0 = 1.5
	
	real :: Rgas=8.3144     ! the universal gas constant J/∫K/mol
	real :: T_leaf
	real :: Gammast, Gamma,gs,Kct, Kot,A,B,C,Vcmaxt,Jmaxt,discriminante,J1,J2,Jk
	real :: IradEf,An,Ci
	real :: X,Gma,Bta,b2,b1,b0,bx,ciquad1,ciquad2,Aquad1,Aquad2
	
	T_leaf = airT
	Ci = Ca*0.8
	
	
	! Bernacchi et al. (2001) temperature dependence of Farquhar et al. 1980 Eq. 38
	Gammast = Gammast25 * exp((EaGammast / (Rgas * 298.1)) * (1 - 298.1 / (T_leaf + 273.1)))
	Kct = Kc * exp((EaKc / (Rgas * 298.1)) * (1 - 298.1 / (T_leaf + 273.1)))
	Kot = Ko * exp((EaKo / (Rgas * 298.1)) * (1 - 298.1 / (T_leaf + 273.1)))
	
	! Vcmax temperature dependence
	A = exp((EaVcmax * (T_leaf + 273.1 - 298.1)) / (298.1 * Rgas * (T_leaf + 273.1)))
	B = 1 + exp((SJmax * 298.1 - EdVcmax) / (Rgas * 298.1))
	C = 1 + exp((SJmax * (273.1 + T_leaf) - EdVcmax) / (Rgas * (273.1 + T_leaf)))
	
	Vcmaxt = Vcmax * A * B / C
	
	! Jmax temperature dependence
	! De Pury and Farquhar et al. 1997   Plant, Cell and Env. 20:537-557
	A = exp((EaJmax * (T_leaf + 273.1 - 298.1)) / (298.1 * Rgas * (T_leaf + 273.1)))
	B = 1 + exp((SJmax * 298.1 - EdJmax) / (Rgas * 298.1))
	C = 1 + exp((SJmax * (273.1 + T_leaf) - EdJmax) / (Rgas * (273.1 + T_leaf)))
	
	! calculate Jmax = f(Vmax), reference:
	! Wullschleger, S.D., 1993.  Biochemical limitations to carbon assimilation
	! in C3 plants - A retrospective analysis of the A/Ci curves from
	! 109 species. Journal of Experimental Botany, 44:907-920.	
	Jmaxt = 2.1*Vcmax * A * B / C
	
	! temperature response of mitochondrial respiration
	Rdt = Rd * (VQ10**((T_leaf - 25) / 10))
	
	! calculate J = f(Jmax, ppfd), reference:
	! de Pury and Farquhar 1997 Plant Cell and Env.
	IradEf=(ppfdACT*pabs)/ppe
	discriminante =sqrt((IradEf  + Jmaxt)**2 - (4 * tetaph * IradEf * Jmaxt))
	J1 = ((IradEf + Jmaxt) + discriminante) / (2 * tetaph)
	J2 = ((IradEf + Jmaxt) - discriminante) / (2 * tetaph)
	jk=min(j1,j2)

	! if there is no activity
	If((Vcmaxt.le.0).or.(Jmaxt.le.0)) Then ! max rate of carboxilation or max rate of electron transport = 0
		An = 0
		gs = gs0    ! rate of photosynthesis = 0, stomatal conductance = gs0
		PhotoSynth = An
		Return
	End If
	
	! gamma is the point of compensation for co2 - the c02 concentration below which results in a
	! negative net photosynthesis
	! Farquhar et al. 1980 Eq. 39
	Gamma = (Gammast + (Kct * (1 + (Oi / Kot))) * Rdt / Vcmaxt) / (1 - (Rdt / Vcmaxt))
	If((Gamma.lt.Gammast).Or.(Gamma.gt.Ca)) Gamma = Gammast

	! Analytical solution for ci. This is the ci which satisfies supply and demand
	! functions simultaneously
	
	!     calculate X using Luening Ball Berry model (should scale for soil moisture?)
	      X = g1/((Ca - Gamma)*(1+VPD/gsD0))
	      	
	!     calculate solution for ci when Rubisco activity limits A
	      Gma = VcmaxT  
	      Bta = KcT*(1.+ Oi/KoT)

	!	calculate coefficients for quadratic equation for ci
		b2 = gs0+X*(Gma-Rdt)
		b1 = (1.-Ca*X)*(Gma-Rdt)+gs0*(Bta-Ca)-X*(Gma*gammast+Bta*Rdt)
		b0 = -(1.-Ca*X)*(Gma*gammast+Bta*Rdt)-gs0*Bta*Ca
		bx=b1*b1-4.*b2*b0
	
	! calculate larger root of quadratic
		if(bx.gt.0.0)ciquad1 = (-b1+sqrt(bx))/(2.*b2)
		
	      	if((ciquad1.lt.0).or.(bx.lt.0.))then
			Aquad1 = 0.0
			ciquad1 = Ci
		else
			Aquad1 = Gma*(ciquad1-gamma)/(ciquad1+Bta)
		endif

	!     calculate +ve root for ci when RuBP regeneration limits A
	      Gma = Jk/4.
	      Bta = 2.*gamma

	!     	calculate coefficients for quadratic equation for ci		
		b2 = gs0+X*(Gma-Rdt)
		b1 = (1.-Ca*X)*(Gma-Rdt)+gs0*(Bta-Ca)-X*(Gma*gammast+Bta*Rdt)
		b0 = -(1.-Ca*X)*(Gma*gammast+Bta*Rdt)-gs0*Bta*Ca
		bx=b1*b1-4.*b2*b0
	
	! 	calculate larger root of quadratic
		if(bx.gt.0.0) ciquad2 = (-b1+sqrt(bx))/(2.*b2)
		
	      	if((ciquad2.lt.0).or.(bx.lt.0.))then
			Aquad2 = 0.0
			ciquad2 = Ci
		else
			Aquad2 = Gma*(ciquad2-gamma)/(ciquad2+Bta)
		endif

	!     choose smaller of Ac, Aq
		An = (min(Aquad1,Aquad2) - Rdt)
	
 	! 	calculate new values for gsc, cs (BBL model) Leuning R. 1995, Plant,Cell and Environ. 18:339-355 and L. et al. 18:1183-1200
		If ((An).gt.0) Then
			gs= gs0 + ((g1 * (An) * 1.6) / ((Ca - Gamma) * (1 + (VPD / gsD0))))
		Else
			gs = gs0
		End If
		
		! Gc = (gs*VPD)*18.01528*10/1000
        Gc = gs * 1.6 * VPD * 1000 / 101300 * 18 / 1000
        ! This Gc is in g H2O m-2 s-1


        ! gs is in mol/m^2 leaf area per day
		! there are 18.01528g water per mol
		! multiply by the density of water 1cm^3/g
		! multipy by 100cm and divide by 10^6cm^3
        ! Gc is now in mm m-2 = kg m-2

		If ((gs - gs0).gt.0) Then
			Ci = Ca - ((((An) * 1.6) / (gs - gs0)))
			If (Ci.lt.Gamma) Then
			! Recalculate Anph to adjust to lowest conductance.
				An = (((Ca - Gamma) * (gs - gs0)) / 1.6) + Rdt
				Ci = Gamma
			End If
		Else
			Ci = Gamma
			An =(((Ca - Gamma) * (gs - gs0)) / 1.6)
		End If
 
		PhotoSynth = An

	return

END
