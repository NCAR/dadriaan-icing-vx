;
; File: match_ij_data.ncl
;
; Author: D. Adriaansen
;
; Date: 23 May 2019
;
; Purpose: Find grid cells in a neighborhood matched to a point observation (lat/lon)
;
; Notes:
;___________________________________________________________________________

; Load NCL libraries
load "/home/dadriaan/cvs/pat/dadriaan/scripts/ncl/load_ncl_libs.ncl"

;##################################### User Config ################################;

; Path to netCDF file with grid information in it
nc = "/home/dadriaan/projects/sae2019/data/grids/hrrr/20190201/g_000000/f_00010800.nc"
;nc = "/home/dadriaan/projects/sae2019/data/grids/rap/20190201/g_000000/f_00010800.nc"

; Path to point data with lat/lon info
ptobs = "/home/dadriaan/projects/sae2019/data/pireps/PIREPSpos.csv"
;ptobs = "/home/dadriaan/projects/sae2019/data/pireps/PIREPSposneg.csv"

; Path to output file
;out = "/home/dadriaan/projects/sae2019/data/match/hrrr/PIREPShrrr.out"
out = "/home/dadriaan/projects/sae2019/data/match/rap/PIREPSrap.out"

; How many points on the side of the square in our neighborhood?
; For HRRR, 8 = 8x8 = 64 pts = 24 km^2 area
; For RAP, 2 = 2x2 = 4 pts = 26 km^2 area
;nside = 8
nside = 2

; Offset, to start at the largest bounding square
offset = (nside/2)-1

; How many total points in our neighborhood
npts = nside*nside

; What forecast lead time (seconds) do we want to use for matching model data to obs?
nsecfcst = 10800
;nsecfcst = 21600

; What time delta do we want for matching model data to obs (seconds)?
matchdt = 1800

; Do we want to match before, after, or both?
;matchwindow = "both"
matchwindow = "after"
;matchwindow = "before"

; Match minutes
; For FIP, set to 0. For CIP, use minutes after the hour that CIP was run.
matchmins = 8

; What mode? Options are "CIP" or "FIP"
mode = "CIP"

;##################################################################################;

print("")
print("WARNING! This code expects a model grid with the lower left corner as point (0,0)")
print("")

; Clean up the old file if present
if fileexists(out) then
  system("rm -rfv "+out)
end if

; Open up the point data file
odata = asciiread(ptobs,-1,"string")

; Read data- customize to input file
ut = toint(str_get_field(odata,1,","))
olat = tofloat(str_get_field(odata,2,","))
olon = tofloat(str_get_field(odata,3,","))
rawObs = str_get_field(odata,35,",")
flvl = toint(str_get_field(odata,4,","))
temp = tofloat(str_get_field(odata,13,","))
ibase1 = toint(str_get_field(odata,16,","))
itop1 = toint(str_get_field(odata,17,","))
iint1 = toint(str_get_field(odata,18,","))
ityp1 = toint(str_get_field(odata,19,","))
ibase2 = toint(str_get_field(odata,20,","))
itop2 = toint(str_get_field(odata,21,","))
iint2 = toint(str_get_field(odata,22,","))
ityp2 = toint(str_get_field(odata,23,","))
acft = str_get_field(odata,34,",")
  
; Open the netCDF file with the grid information
fin = addfile(nc,"r")
lat2d = fin->lat0
lon2d = fin->lon0

; Find number of points in grid
gDims = dimsizes(lat2d)
numY = gDims(0)-1
numX = gDims(1)-1

; Counter for obs
ocnt = 1

; Counter for processed obs
idcnt = 0

; Loop over each observation in the file
do n=0,dimsizes(ut)-1,1
;do n=0,10,1

  ;print("")
  ;print("EXAMINING TIME: "+ut(n))
  ;print("")
 
  print("")
  print("PROCESSING OB: "+rawObs(n)+" ("+sprinti("%.5d",toint(ocnt))+"/"+sprinti("%.5d",toint(dimsizes(rawObs)))+")")
  print("")

  ; Current time
  curr = ut(n)

  ; Convert UNIX observation time to time elements (YMDHMS)
  curr@units = "seconds since 1970-1-1 00:00:00"
  obtime = cd_calendar(curr,0)
  oyr = obtime(:,0)
  omo = obtime(:,1)
  ody = obtime(:,2)
  ohr = obtime(:,3)
  omn = obtime(:,4)
  osc = obtime(:,5)

  ; Create a lower time equal to HR:00:00 of the observation. This is essentially floor(obtime), or the top of the hour
  ; containing the observation. tlower always <= obtime.
  opt = 0
  tlower = cd_inv_calendar(toint(oyr),toint(omo),toint(ody),toint(ohr),toint(matchmins),0,"seconds since 1970-1-1 00:00:00",opt)
  if tlower.gt.curr then
    tlower = tlower-3600
  end if

  ; Create an upper time equal to HR+1:00:00 of the observation. This will be the upper bound of the observation time.
  ; tupper always >= obtime.
  tupper = tlower+3600

  ; Figure out the tmatch, based on the matching method requested
  if matchwindow.eq."both" then

    ; Find the closest time. If the delta between the observation and the lower time is < 30 minutes, then it will
    ; match to the previous top of the hour. Otherwise (delta >= 30 minutes), it will jump to the next hour.
    if (curr-tlower).lt.matchdt then
      tmatch = tlower
    else
      tmatch = tupper
    end if
  end if
  if matchwindow.eq."after" then
    
    ; Find the closest time. If the delta between the observation and the lower time is < matchdt seconds, then it will
    ; match to the tlower. Othwerwise (delta >= matchdt seconds), it will skip this PIREP and not use it.
    if (curr-tlower).lt.matchdt then
      tmatch = tlower
    else
      print("")
      print("SKIPPING PIREP. OUTSIDE TIME BOUNDS.")
      continue
    end if
  end if
  if matchwindow.eq."before" then
  
    ; Find the closest time. If the delta between the observation and the upper time is < matchdt seconds, then it will
    ; match to the tupper. Otherwise (delta >= matchdt seconds), it will skip this PIREP and not use it.
    if (tupper-curr).lt.matchdt then
      tmatch = tupper
    else
      print("")
      print("SKIPPING PIREP. OUTSIDE THE TIME BOUNDS.")
      continue
    end if
  end if
  ;print(""+matchwindow)
  ;print(""+matchdt)
  ;print(""+matchmins)
  ;cStr = prettyPrint(curr,"ymdhms")
  ;print(""+cStr)
  ;mStr = prettyPrint(toint(tmatch),"ymdhms")
  ;print(""+mStr)
  ;continue 

  ; Back out the init time/date for the best match
  if mode.eq."FIP" then
    init = tmatch-(nsecfcst)
  else
    init = tmatch
  end if
  init@units = "seconds since 1970-1-1 00:00:00"
  di = cd_calendar(init,0)
  iyr = di(:,0)
  imo = di(:,1)
  idy = di(:,2)
  ihr = di(:,3)
  istring = sprinti("%.4i",toint(iyr))+sprinti("%.2i",toint(imo))+sprinti("%.2i",toint(idy))
  if mode.eq."FIP" then
    gstring = "g_"+sprinti("%.2i",toint(ihr))+"0000"
    fname = "f_"+sprinti("%.8i",toint(nsecfcst))+".nc"
    file_string = istring+"/"+gstring+"/"+fname
  else
    fname = sprinti("%.2i",toint(ihr))+sprinti("%.2i",toint(matchmins))+"00.nc"
    file_string = istring+"/"+istring+"_"+fname
  end if

  ; Construct a string identifying the path to the correct model file for these times
  ; TODO: Only do this when we are finding exact files. We can just write out the istring,gstring,fname and match to
  ; parentDir later in another script

  ; Locate the closest grid cell to this observation
  ij = getind_latlon2d(lat2d,lon2d,olat(n),olon(n))
  j = ij(0,0)
  i = ij(0,1)
  homeI = i
  homeJ = j
  
  print("")
  print("NEAREST POINT TO OB:")
  print("i = "+i)
  print("j = "+j)

  ; Determine the neighborhood of grid cells for this point. The initial i,j is the closest point to the observation so to
  ; determine which 4 points "encapsulate" the observation we need to check the great
  ; circle distances between the observation and the following four points:
  ; (j+1,i+1) --> UR
  ; (j+1,i-1) --> UL
  ; (j-1,i+1) --> LR
  ; (j-1,i-1) --> LL
  ; Whichever of those four distances is closest represents the opposite corner of the four "encapsulating" points and
  ; we can construct the first ring of four points around the observation.
  ; NOTE: The starting i,j will be the northwest (upper left) corner of the four bounding points. Looping occurs XXXXXX
  URi = i+1
  URj = j+1
  if URi.gt.numX .or. URi.lt.0 .or. URj.gt.numY .or. URj.lt.0 then
    print("")
    print("FATAL! UR corner off grid. Unable to process OB "+rawObs(n)+".")
    print("")
    ocnt = ocnt + 1
    continue
  else
    UR = gc_latlon(olat(n),olon(n),lat2d(URj,URi),lon2d(URj,URi),10,4)
  end if
  LRi = i+1
  LRj = j-1
  if LRi.gt.numX .or. LRi.lt.0 .or. LRj.gt.numY .or. LRj.lt.0 then
    print("")
    print("FATAL! LR corner off grid. Unable to process OB "+rawObs(n)+".")
    print("")
    ocnt = ocnt + 1
    continue
  else
    LR = gc_latlon(olat(n),olon(n),lat2d(LRj,LRi),lon2d(LRj,LRi),10,4)
  end if
  ULi = i-1
  ULj = j+1
  if ULi.gt.numX .or. ULi.lt.0 .or. ULj.gt.numY .or. ULj.lt.0 then
    print("")
    print("FATAL! UL corner off grid. Unable to process OB "+rawObs(n)+".")
    print("")
    ocnt = ocnt + 1 
    continue
  else
    UL = gc_latlon(olat(n),olon(n),lat2d(ULj,ULi),lon2d(ULj,ULi),10,4)
  end if
  LLi = i-1
  LLj = j-1
  if LLi.gt.numX .or. LLi.lt.0 .or. LLj.gt.numY .or. LLj.lt.0 then
    print("")
    print("FATAL! LL corner off grid. Unable to process OB "+rawObs(n)+".")
    print("")
    ocnt = ocnt + 1
    continue
  else
    LL = gc_latlon(olat(n),olon(n),lat2d(LLj,LLi),lon2d(LLj,LLi),10,4)
  end if
  dChk = (/UL,LL,UR,LR/)
  ;box = ind(min((/UL,LL,UR,LR/))) ; BUG!
  box = ind(dChk.eq.min(dChk))
  print("")
  print("DISTANCES: ")
  print("UL ("+ULi+","+ULj+") dist = "+dChk(0))
  print("LL ("+LLi+","+LLj+") dist = "+dChk(1))
  print("UR ("+URi+","+URj+") dist = "+dChk(2))
  print("LR ("+LRi+","+LRj+") dist = "+dChk(3))
  print("")
  corner = ""
  ; UPPER LEFT
  if box.eq.0 then
    print("FOUND UL ENCAPSULATION")
    print("NEAREST IS LR CORNER OF ENCAPSULATION")
    ;si = (i-1)-offset
    ;sj = (j+1)+offset
    minI = (i-1)-offset ; Left boundary
    minJ = j-offset ; Bottom boundary
    maxI = (minI+nside)-1 ; Right boundary
    maxJ = (minJ+nside)-1 ; Top boundary
    ULi = i-1
    ULj = j+1
    LLi = i-1
    LLj = j
    URi = i
    URj = j+1
    LRi = i
    LRj = j
    print("ULe = "+ULi+" "+ULj)
    print("URe = "+URi+" "+URj)
    print("LRe = "+LRi+" "+LRj)
    print("LLe = "+LLi+" "+LLj)
    boxfound = True
    corner = "LR"
  end if
  ; LOWER LEFT
  if box.eq.1 then
    print("FOUND LL ENCAPSULATION")
    print("NEAREST IS UR CORNER OF ENCAPSULATION")
    ;si = (i-1)-offset
    ;sj = j+offset
    minI = (i-1)-offset ; Left boundary
    minJ = (j-1)-offset ; Bottom boundary
    maxI = (minI+nside)-1 ; Right boundary
    maxJ = (minJ+nside)-1 ; Top boundary
    ULi = i-1
    ULj = j
    LLi = i-1
    LLj = j-1
    URi = i
    URj = j
    LRi = i
    LRj = j-1
    print("ULe = "+ULi+" "+ULj)
    print("URe = "+URi+" "+URj)
    print("LRe = "+LRi+" "+LRj)
    print("LLe = "+LLi+" "+LLj)
    boxfound = True
    corner = "UR"
  end if
  ; UPPER RIGHT
  if box.eq.2 then
    print("FOUND UR ENCAPSULATION")
    print("NEAREST IS LL CORNER OF ENCAPSULATION")
    ;si = (i+1)-offset
    ;sj = (j+1)+offset
    minI = i-offset ; Left boundary
    minJ = (j-1)-offset ; Bottom boundary
    maxI = (minI+nside)-1 ; Right boundary
    maxJ = (minJ+nside)-1 ; Top boundary
    ULi = i
    ULj = j+1
    LLi = i
    LLj = j
    URi = i+1
    URj = j+1
    LRi = i+1
    LRj = j
    print("ULe = "+ULi+" "+ULj)
    print("URe = "+URi+" "+URj)
    print("LRe = "+LRi+" "+LRj)
    print("LLe = "+LLi+" "+LLj)
    boxfound = True
    corner = "LL"
  end if
  ; LOWER RIGHT
  if box.eq.3 then
    print("FOUND LR ENCAPSULATION")
    print("NEAREST IS UL CORNER OF ENCAPSULATION")
    ;si = (i-offset)
    ;sj = (j+offset)
    minI = i-offset ; Left boundary
    minJ = j-offset ; Botton boundary
    maxI = (minI+nside)-1 ; Right boundary
    maxJ = (minJ+nside)-1 ; Top boundary
    ULi = i
    ULj = j
    LLi = i
    LLj = j-1
    URi = i+1
    URj = j
    LRi = i+1
    LRj = j-1
    print("ULe = "+ULi+" "+ULj)
    print("URe = "+URi+" "+URj)
    print("LRe = "+LRi+" "+LRj)
    print("LLe = "+LLi+" "+LLj)
    boxfound = True
    corner = "UL"
  end if
  
  ; Sanity check
  if .not.boxfound then
    print("")
    print("FATAL! BOUNDING BOX NOT FOUND FOR OB: "+ut(n))
    print("")
    exit()
  end if


  ; Make sure the minI, minJ, maxI, maxJ are "in bounds"
  if minI.lt.0 .or. minI.gt.numX .or. maxI.gt.numX .or. maxI.lt.0 .or. minJ.lt.0 .or. minJ.gt.numY .or. maxJ.lt.0 .or. maxJ.gt.numY then
    print("")
    print("FATAL! Neighborhood extends outside grid. Unable to process OB "+rawObs(n)+".")
    print("")
    ocnt = ocnt + 1
    continue
  end if

  ; Print out info to see if it makes sense
  print("")
  print("BOUNDARIES OF ALL POINTS IN NEIGHBORHOOD:")
  print("minI = "+minI)
  print("minJ = "+minJ)
  print("maxI = "+maxI)
  print("maxJ = "+maxJ)
  print("################################################")

  ; Total number of points using minI/maxI minJ/maxJ
  numI = dimsizes(ispan(minI,maxI,1))
  numJ = dimsizes(ispan(minJ,maxJ,1))
  totPts = numI*numJ

  ; Write out some data to a file
  ; ID,unix_time (obs),lat,lon,flvl,tmp,ibase1,itop1,iint1,ityp1,ibase2,itop2,iint2,ityp2,acft,rawObs,unix_time (valid),forecast,minI,maxI,minJ,maxJ,npts,file_string
  format = "%d,%d,%f,%f,%d,%f,%d,%d,%d,%d,%d,%d,%d,%d,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%s,%d,%s"
  alist = [/idcnt,ut(n),olat(n),olon(n),flvl(n),temp(n),ibase1(n),itop1(n),iint1(n),ityp1(n),ibase2(n),itop2(n),iint2(n),ityp2(n),acft(n),rawObs(n),toint(tmatch),nsecfcst,minI,maxI,minJ,maxJ,homeI,homeJ,corner,totPts,file_string/]
  write_table(out,"a",alist,format)

  ; Advance counters
  ocnt = ocnt + 1
  idcnt = idcnt + 1

end do
