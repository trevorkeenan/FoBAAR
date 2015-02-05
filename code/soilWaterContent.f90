REAL FUNCTION soilWaterContent(soilWater,precipDaily,ETdaily,swhc,drainageP,drainage,runoff) 

implicit none
  	real, intent(in) :: soilWater,swhc,precipDaily,ETdaily,drainageP
	real, intent(out) :: drainage,runoff
	 
	! calculate the soil water content as
	! SWC(i+1)=SWC(i)+precip-ET-drainage-runoff
	! runoff occurs when SWC=soil water holding capacity (swhc)
	! drainage is a constant fraction of SWC
	soilWaterContent = soilWater+precipDaily-ETdaily

        if (soilWaterContent.lt.0)then
                soilWaterContent=0
        endif
        
        ! calculate runoff
        if (soilWaterContent.gt.swhc)then
                runoff=soilWaterContent-swhc
        else
                runoff=0
        endif
        soilWaterContent=soilWaterContent-runoff
        
        ! calculate drainage
        drainage=drainageP*soilWaterContent
        soilWaterContent=soilWaterContent-drainage
        
	return

END
