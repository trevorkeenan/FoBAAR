REAL FUNCTION TotalErrorV4bayes(err,RealityErr,rep,numConstraints,constraints,numSoilPools,toterrAll,countD,nparams,boundsP,P)
! this function will accept the error matrix and calulate total model data mismatch
! for the constraints indicated in the constraints matrix.

! the cost function includes the parameter values and priors

implicit none
        integer, intent(in) :: numConstraints,numSoilPools,rep,nparams
	
  	real, intent(in) :: countD(numConstraints)
  	real, intent(in) :: constraints(numConstraints+1)
	real, intent(in) :: RealityErr
	real,intent(in) :: boundsP(nparams,2),P(nparams)

        double PRECISION, intent(out) :: toterrAll
        real, intent(out) :: err(numConstraints,2)
	real :: tmp2 
        integer :: i
        logical, dimension(numConstraints+1) :: mask
               

       
	do i =1,numConstraints
				if(i.ne.9)then	! 9 is long-term increment data - no average
				! get average error for each data stream
				       if(err(i,1).gt.0)err(i,1)=err(i,1)/countD(i)  
				endif
	enddo
			
	do i=1,numConstraints 
			! check for NaN errors
				if(err(i,1).lt.10E16)then
				else !it must be -NaN
					err(i,1)=10E16
				endif
	enddo
		
        ! set the number of constraints to be used
	if (numSoilPools.eq.1)then
		tmp2=20
	else if (numSoilPools.eq.2)then
		tmp2=22
	else
		tmp2=24
	endif
	
	mask=constraints.eq.1
	tmp2=count(mask,1)
	
	if(mask(9))tmp2=tmp2-1
	if(mask(18))tmp2=tmp2-1
	if(mask(28).and.(numSoilPools.eq.1))tmp2=tmp2-1
	
	! calculate the total error considering all constraints
	! for observations
	 TotalErrorV4bayes=0
		        do i=1,numConstraints
		        	! loop through error matrix. If we're using a data stream as a constraint, then add to total error value
		        	TotalErrorV4bayes = TotalErrorV4bayes+err(i,1)*constraints(i)
		        	
			enddo
			TotalErrorV4bayes = TotalErrorV4bayes+RealityErr*constraints(i)
			
			TotalErrorV4bayes =TotalErrorV4bayes/tmp2
                	! for parameters
		        toterrAll=0
		        do i=1,nparams
		        	! loop through error matrix. If we're using a data stream as a constraint, then add to total error value
		        	tmp2=(P(i)- ((boundsP(i,1)+boundsP(i,2))/2))**2/(((boundsP(i,1)+boundsP(i,2)))**2)
		        	toterrAll = toterrAll+tmp2
		        	
			enddo
			toterrAll=toterrAll/nparams
		        
		        TotalErrorV4bayes=1000*(TotalErrorV4bayes+toterrAll)/2
        	      
			if (TotalErrorV4bayes.le.0) TotalErrorV4bayes = 10E10
			if(TotalErrorV4bayes.lt.10E10)then
			else !it must be -NaN
				TotalErrorV4bayes=10E10
			endif

	return

END

