C     Diagram 1
      DATA MAPCONFIG(1)/1/
      DATA TSTRATEGY(1)/1/
      DATA (IFOREST(I,-1,1),I=1,2)/2,5/
      DATA TPRID(-1,1)/24/
      DATA (SPROP(I,-1,1),I=1,8)/0,0,0,0,0,0,0,0/
      DATA (IFOREST(I,-2,1),I=1,2)/-1,3/
      DATA TPRID(-2,1)/24/
      DATA (SPROP(I,-2,1),I=1,8)/0,0,0,0,0,0,0,0/
      DATA (IFOREST(I,-3,1),I=1,2)/-2,4/
      DATA TPRID(-3,1)/2/
      DATA (SPROP(I,-3,1),I=1,8)/0,0,0,0,0,0,0,0/
C     Diagram 2
      DATA MAPCONFIG(2)/2/
      DATA TSTRATEGY(2)/2/
      DATA (IFOREST(I,-1,2),I=1,2)/1,4/
      DATA TPRID(-1,2)/24/
      DATA (SPROP(I,-1,2),I=1,8)/0,0,0,0,0,0,0,0/
      DATA (IFOREST(I,-2,2),I=1,2)/-1,5/
      DATA TPRID(-2,2)/4/
      DATA (SPROP(I,-2,2),I=1,8)/0,0,0,0,0,0,0,0/
      DATA (IFOREST(I,-3,2),I=1,2)/-2,3/
      DATA TPRID(-3,2)/4/
      DATA (SPROP(I,-3,2),I=1,8)/0,0,0,0,0,0,0,0/
C     Diagram 3
      DATA MAPCONFIG(3)/3/
      DATA TSTRATEGY(3)/2/
      DATA (IFOREST(I,-1,3),I=1,2)/5,4/
      DATA (SPROP(I,-1,3),I=1,8)/0,0,0,0,-24,-24,-24,-24/
      DATA TPRID(-1,3)/0/
      DATA (IFOREST(I,-2,3),I=1,2)/-1,3/
      DATA (SPROP(I,-2,3),I=1,8)/0,0,0,0,-24,-24,-24,-24/
      DATA TPRID(-2,3)/0/
C     Diagram 4
      DATA MAPCONFIG(4)/4/
      DATA TSTRATEGY(4)/2/
      DATA (IFOREST(I,-1,4),I=1,2)/5,4/
      DATA (SPROP(I,-1,4),I=1,8)/0,0,0,0,-24,-24,-24,-24/
      DATA TPRID(-1,4)/0/
      DATA (IFOREST(I,-2,4),I=1,2)/1,-1/
      DATA TPRID(-2,4)/4/
      DATA (SPROP(I,-2,4),I=1,8)/0,0,0,0,0,0,0,0/
      DATA (IFOREST(I,-3,4),I=1,2)/-2,3/
      DATA TPRID(-3,4)/4/
      DATA (SPROP(I,-3,4),I=1,8)/0,0,0,0,0,0,0,0/
C     Number of configs
      DATA MAPCONFIG(0)/4/
