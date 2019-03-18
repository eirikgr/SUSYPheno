
	if( sqrtm .lt. 11.832159566199232D0 )
     &    Warning('Extrapolating tHm in MHiggs')
	if( TBeff .lt. 2.5D0 )
     &    Warning('Extrapolating tHm in TBeff')

	if( sqrtm.lt.24.49489742783178D0 ) then

        tHm = exp(-11387.267296834201D0 + 
     &    TBeff*(9.265665516009038D0 + 
     &       TBeff*(0.1287355144336406D0 + 
     &          TBeff*(-0.015546423693133814D0 + 
     &             TBeff*(0.0006312217906869511D0 + 
     &                TBeff*
     &                 (-0.000016726287404224414D0 + 
     &                   TBeff*
     &                    (2.8372875311555307D-7 + 
     &                      TBeff*
     &                       (-2.988567462573199D-9 + 
     &                       TBeff*
     &                       (1.785649085682255D-11 - 
     &                       4.6290905426189013D-14*TBeff))))))))+
     &      sqrtm*(5744.576094117645D0 + 
     &       TBeff*(-4.60454516731925D0 + 
     &          TBeff*(0.03854216548588106D0 + 
     &             TBeff*(0.0001726839749845425D0 + 
     &                TBeff*
     &                 (-1.4912612603549435D-6 + 
     &                   TBeff*
     &                    (2.9655521188258198D-8 + 
     &                      TBeff*
     &                       (-2.023954501864769D-10 + 
     &                       TBeff*
     &                       (-1.455472486602509D-12 + 
     &                       1.560517412823355D-14*TBeff))))))) + 
     &       sqrtm*(-1277.123400324312D0 + 
     &          TBeff*(0.8329154446883827D0 + 
     &             TBeff*(-0.0070257233781612135D0 + 
     &                TBeff*
     &                 (-0.00001975070351237301D0 + 
     &                   TBeff*
     &                    (2.7073110842621367D-8 + 
     &                      TBeff*
     &                       (-9.124491291178447D-10 + 
     &                       TBeff*
     &                       (1.5861409286707272D-11 - 
     &                       6.087876884093844D-14*TBeff)))))) + 
     &          sqrtm*(164.43475130397033D0 + 
     &             TBeff*(-0.08387125461798306D0 + 
     &                TBeff*
     &                 (0.0006949376628685907D0 + 
     &                   TBeff*
     &                    (1.4308196335324948D-6 + 
     &                      TBeff*
     &                       (6.093309988661061D-10 + 
     &                       TBeff*
     &                       (-2.3389965333257225D-11 - 
     &                       2.313795556953886D-14*TBeff))))) + 
     &             sqrtm*(-13.516640053996065D0 + 
     &                TBeff*
     &                 (0.00512334446398382D0 + 
     &                   TBeff*
     &                    (-0.00004066684695644574D0 + 
     &                      TBeff*
     &                       (-5.973748442932111D-8 + 
     &                       TBeff*
     &                       (4.4096563593978973D-11 + 
     &                       3.8077141934363574D-13*TBeff)))) + 
     &                sqrtm*
     &                 (0.7358036721373993D0 + 
     &                   TBeff*
     &                    (-0.0001932037560623426D0 + 
     &                      TBeff*
     &                       (1.409055910645403D-6 + 
     &                       TBeff*
     &                       (1.21988146872467D-9 - 
     &                       1.2298053618850378D-12*TBeff))) + 
     &                   sqrtm*
     &                    (-0.026532567690852046D0 + 
     &                      TBeff*
     &                       (4.3478273620635665D-6 + 
     &                       TBeff*
     &                       (-2.6696499284740174D-8 - 
     &                       9.307367977481643D-12*TBeff)) + 
     &                      sqrtm*
     &                       (0.0006112624628154316D0 + 
     &                       TBeff*
     &                       (-5.243572309038002D-8 + 
     &                       2.1293710846748402D-10*TBeff) + 
     &                       sqrtm*
     &                       (-8.165927211817924D-6 + 
     &                       4.820628899214805D-8*sqrtm + 
     &                       2.501321855254812D-10*TBeff)))))))))

	else

        tHm = exp(-203034.88582426897D0 + 
     &    TBeff*(686.0532744039647D0 + 
     &       TBeff*(3.7196424090713474D0 + 
     &          TBeff*(-0.04688105344911253D0 + 
     &             TBeff*(0.0004855852633036865D0 + 
     &                TBeff*
     &                 (-0.00001204984954853224D0 + 
     &                   TBeff*
     &                    (1.7919495230365673D-7 + 
     &                      TBeff*
     &                       (-1.4696041925467152D-9 + 
     &                       5.129198947532139D-12*TBeff))))))) + 
     &    sqrtm*(58173.88886103328D0 + 
     &       TBeff*(-176.38493236937714D0 + 
     &          TBeff*(-0.6333160355726964D0 + 
     &             TBeff*(0.0061495420274382175D0 + 
     &                TBeff*
     &                 (3.0339112293440392D-6 + 
     &                   TBeff*
     &                    (-3.6480596930798124D-8 + 
     &                      TBeff*
     &                       (2.6208613247713043D-10 - 
     &                       9.760546057728793D-13*TBeff)))))) + 
     &       sqrtm*(-7286.750963373372D0 + 
     &          TBeff*(19.22718260097787D0 + 
     &             TBeff*(0.04541833050510523D0 + 
     &                TBeff*
     &                 (-0.00044718405941661317D0 + 
     &                   TBeff*
     &                    (-7.697139358358937D-8 + 
     &                      TBeff*
     &                       (4.467457328790491D-10 - 
     &                       8.66175158369903D-13*TBeff))))) + 
     &          sqrtm*(521.2238218229023D0 + 
     &             TBeff*(-1.154652085869666D0 + 
     &                TBeff*
     &                 (-0.001618691523200203D0 + 
     &                   TBeff*
     &                    (0.000016130243283800877D0 + 
     &                      TBeff*
     &                       (1.2187129994619349D-9 - 
     &                       3.0475299730495255D-12*TBeff)))) + 
     &             sqrtm*(-23.28895070443708D0 + 
     &                TBeff*
     &                 (0.04125630544489795D0 + 
     &                   TBeff*
     &                    (0.0000285406245812878D0 + 
     &                      TBeff*
     &                       (-2.899278663709154D-7 - 
     &                       8.701820888484632D-12*TBeff))) + 
     &                sqrtm*
     &                 (0.665643604596355D0 + 
     &                   TBeff*
     &                    (-0.0008769769351036599D0 + 
     &                      TBeff*
     &                       (-1.9398395068437547D-7 + 
     &                       2.0801612766854388D-9*TBeff)) + 
     &                   sqrtm*
     &                    (-0.011885854917238039D0 + 
     &                      TBeff*
     &                       (0.000010266547785091555D0 - 
     &                       1.224909347489249D-10*TBeff) + 
     &                      sqrtm*
     &                       (0.000121234873254401D0 - 
     &                       5.408508995853367D-7*sqrtm - 
     &                       5.1043282696247345D-8*TBeff))))))))

	endif

#ifdef DETAILED_DEBUG
	DPROD 'tHm =', tHm ENDL
#endif

