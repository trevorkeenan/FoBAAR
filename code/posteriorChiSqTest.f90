REAL FUNCTION posteriorChiSqTest(err,forest,rep,decidflag,RealityErr,numConstraints,constraints)
! This is the betaV4 posterior chi-sq test
! changes: 
! - predefined cut-off array
! - only tests for constraints listed in constraints file

implicit none
  	integer, intent(in) :: numConstraints,decidflag,rep

        real, intent(in) :: RealityErr
  	real, intent(in) :: err(numConstraints,2)
  	real, intent(in) :: constraints(numConstraints+1)
	real :: cutoff(numConstraints)
        character(len=*), intent(in) :: forest
 	integer :: i
        integer :: flag2
	
	cutoff=(/1.06,1.06,1.5,1.5,1.45,&         ! 5
			&1.03,1.0,1.35,100.0,1.05,&   ! 10
			&1.79,2.49,1.4,1.7,1.05,&     ! 15
			&1.2,1.2,1.0,1.06,1.02,&         ! 20
			&1.02,1.0,1.02,1.10,1.1,&     ! 25   
			&1.10,1.00,1.05/)
	
	flag2=1
!	if (decidflag.eq.1)then
		do i=1,numConstraints
			if (constraints(i).eq.1)then
				if (err(i,2).ge.cutoff(i))then
					flag2=0	! is the error is outside acceptable range, flag to reject
				endif
			endif
		enddo
!	endif
	if (RealityErr*constraints(29).gt.0.1)flag2=0	! realityErr is constraint 29

	

posteriorChiSqTest=flag2
return
END


	
	