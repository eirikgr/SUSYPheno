      subroutine dsrddof150
c_______________________________________________________________________
c     table of effective degrees of freedom in the early universe
c
c     t [gev]    g_{\star}^{1/2}   g_{entropy}   for  t_{qcd} = 150 mev
c
c     common:
c       'dsrdcom.h' - included common blocks
c
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c=======================================================================
      implicit none
      include 'dsrdcom.h'
      real*8 table(3,300)
      integer i
      integer iq, ir
c-----------------------------------------------------------------------
      data ((table(iq,ir),iq=1,3),ir=1,99)/
     & 0.1258925d5,10.283486d0,105.749901d0,
     & 0.1188504d5,10.283488d0,105.749901d0,
     & 0.1122019d5,10.283484d0,105.749802d0,
     & 0.1059256d5,10.283484d0,105.749802d0,
     & 0.1000000d5,10.283486d0,105.749802d0,
     & 0.9440606d4,10.283487d0,105.749802d0,
     & 0.8912530d4,10.283480d0,105.749702d0,
     & 0.8413950d4,10.283487d0,105.749702d0,
     & 0.7943284d4,10.283490d0,105.749702d0,
     & 0.7498942d4,10.283487d0,105.749603d0,
     & 0.7079457d4,10.283490d0,105.749603d0,
     & 0.6683455d4,10.283483d0,105.749496d0,
     & 0.6309586d4,10.283490d0,105.749496d0,
     & 0.5956636d4,10.283485d0,105.749397d0,
     & 0.5623426d4,10.283485d0,105.749298d0,
     & 0.5308843d4,10.283485d0,105.749199d0,
     & 0.5011873d4,10.283497d0,105.749199d0,
     & 0.4731512d4,10.283500d0,105.749100d0,
     & 0.4466837d4,10.283503d0,105.749001d0,
     & 0.4216965d4,10.283496d0,105.748802d0,
     & 0.3981071d4,10.283501d0,105.748703d0,
     & 0.3758375d4,10.283501d0,105.748497d0,
     & 0.3548141d4,10.283503d0,105.748299d0,
     & 0.3349663d4,10.283502d0,105.748100d0,
     & 0.3162285d4,10.283512d0,105.747902d0,
     & 0.2985382d4,10.283519d0,105.747704d0,
     & 0.2818396d4,10.283517d0,105.747398d0,
     & 0.2660725d4,10.283524d0,105.747101d0,
     & 0.2511887d4,10.283527d0,105.746696d0,
     & 0.2371374d4,10.283530d0,105.746300d0,
     & 0.2238721d4,10.283529d0,105.745796d0,
     & 0.2113489d4,10.283541d0,105.745300d0,
     & 0.1995262d4,10.283551d0,105.744797d0,
     & 0.1883650d4,10.283556d0,105.744102d0,
     & 0.1778279d4,10.283563d0,105.743401d0,
     & 0.1678808d4,10.283572d0,105.742599d0,
     & 0.1584893d4,10.283583d0,105.741699d0,
     & 0.1496235d4,10.283595d0,105.740700d0,
     & 0.1412538d4,10.283612d0,105.739601d0,
     & 0.1333521d4,10.283621d0,105.738297d0,
     & 0.1258929d4,10.283642d0,105.736900d0,
     & 0.1188502d4,10.283661d0,105.735298d0,
     & 0.1122018d4,10.283676d0,105.733498d0,
     & 0.1059254d4,10.283702d0,105.731499d0,
     & 0.1000000d4,10.283710d0,105.729202d0,
     & 0.9440607d3,10.283741d0,105.726700d0,
     & 0.8912530d3,10.283781d0,105.723801d0,
     & 0.8413969d3,10.283800d0,105.720596d0,
     & 0.7943284d3,10.283820d0,105.717003d0,
     & 0.7498942d3,10.283851d0,105.712997d0,
     & 0.7079456d3,10.283890d0,105.708504d0,
     & 0.6683440d3,10.283950d0,105.703499d0,
     & 0.6309573d3,10.284021d0,105.697800d0,
     & 0.5956623d3,10.284051d0,105.691399d0,
     & 0.5623413d3,10.284140d0,105.684303d0,
     & 0.5308843d3,10.284220d0,105.676201d0,
     & 0.5011873d3,10.284279d0,105.667198d0,
     & 0.4731523d3,10.284350d0,105.657097d0,
     & 0.4466847d3,10.284448d0,105.645798d0,
     & 0.4216965d3,10.284580d0,105.633102d0,
     & 0.3981071d3,10.284698d0,105.618797d0,
     & 0.3758391d3,10.284810d0,105.602798d0,
     & 0.3548134d3,10.284950d0,105.584900d0,
     & 0.3349655d3,10.285129d0,105.564796d0,
     & 0.3162278d3,10.285310d0,105.542198d0,
     & 0.2985382d3,10.285480d0,105.516899d0,
     & 0.2818383d3,10.285710d0,105.488602d0,
     & 0.2660725d3,10.285929d0,105.456802d0,
     & 0.2511887d3,10.286171d0,105.421303d0,
     & 0.2371374d3,10.286439d0,105.381500d0,
     & 0.2238727d3,10.286731d0,105.336998d0,
     & 0.2113489d3,10.287050d0,105.287201d0,
     & 0.1995262d3,10.287380d0,105.231499d0,
     & 0.1883658d3,10.287721d0,105.169296d0,
     & 0.1778279d3,10.288090d0,105.099899d0,
     & 0.1678804d3,10.288471d0,105.022400d0,
     & 0.1584893d3,10.288869d0,104.936096d0,
     & 0.1496236d3,10.289300d0,104.839897d0,
     & 0.1412538d3,10.289671d0,104.732803d0,
     & 0.1333521d3,10.290009d0,104.613899d0,
     & 0.1258925d3,10.290360d0,104.481903d0,
     & 0.1188505d3,10.290609d0,104.335503d0,
     & 0.1122021d3,10.290811d0,104.173599d0,
     & 0.1059254d3,10.290900d0,103.994499d0,
     & 0.1000000d3,10.290780d0,103.796997d0,
     & 0.9440609d2,10.290480d0,103.579597d0,
     & 0.8912511d2,10.289920d0,103.340698d0,
     & 0.8413950d2,10.289018d0,103.078796d0,
     & 0.7943282d2,10.287670d0,102.792397d0,
     & 0.7498942d2,10.285829d0,102.480103d0,
     & 0.7079458d2,10.283350d0,102.140404d0,
     & 0.6683440d2,10.280070d0,101.772202d0,
     & 0.6309573d2,10.275870d0,101.374397d0,
     & 0.5956636d2,10.270591d0,100.946198d0,
     & 0.5623413d2,10.264030d0,100.487000d0,
     & 0.5308845d2,10.255980d0,99.996727d0,
     & 0.5011873d2,10.246230d0,99.475601d0,
     & 0.4731512d2,10.234546d0,98.924393d0,
     & 0.4466836d2,10.220711d0,98.344353d0/
      data ((table(iq,ir),iq=1,3),ir=100,198)/
     & 0.4216965d2,10.204457d0,97.737350d0,
     & 0.3981072d2,10.185567d0,97.105957d0,
     & 0.3758374d2,10.163840d0,96.453339d0,
     & 0.3548133d2,10.139090d0,95.783348d0,
     & 0.3349654d2,10.111179d0,95.100487d0,
     & 0.3162278d2,10.080021d0,94.409843d0,
     & 0.2985389d2,10.045588d0,93.717056d0,
     & 0.2818383d2,10.007990d0,93.028183d0,
     & 0.2660725d2,9.967407d0,92.349503d0,
     & 0.2511886d2,9.924076d0,91.687370d0,
     & 0.2371374d2,9.878415d0,91.048103d0,
     & 0.2238721d2,9.830911d0,90.437637d0,
     & 0.2113489d2,9.782166d0,89.861427d0,
     & 0.1995262d2,9.732853d0,89.324127d0,
     & 0.1883649d2,9.683722d0,88.829483d0,
     & 0.1778279d2,9.635539d0,88.380142d0,
     & 0.1678804d2,9.589060d0,87.977592d0,
     & 0.1584893d2,9.545019d0,87.622047d0,
     & 0.1496235d2,9.504038d0,87.312538d0,
     & 0.1412541d2,9.466631d0,87.046928d0,
     & 0.1333521d2,9.433184d0,86.822151d0,
     & 0.1258925d2,9.403888d0,86.634277d0,
     & 0.1188502d2,9.378805d0,86.478897d0,
     & 0.1122018d2,9.357838d0,86.351212d0,
     & 0.1059254d2,9.340742d0,86.246307d0,
     & 0.1000000d2,9.328511d0,86.159431d0,
     & 0.8912509d1,9.309928d0,86.022141d0,
     & 0.7943287d1,9.300078d0,85.908730d0,
     & 0.7079464d1,9.295743d0,85.796242d0,
     & 0.6309576d1,9.294460d0,85.668533d0,
     & 0.5623413d1,9.294641d0,85.514008d0,
     & 0.5011872d1,9.295348d0,85.323357d0,
     & 0.4466839d1,9.296071d0,85.087929d0,
     & 0.3981072d1,9.296438d0,84.798843d0,
     & 0.3548136d1,9.296024d0,84.446854d0,
     & 0.3162278d1,9.294300d0,84.022621d0,
     & 0.2818383d1,9.290553d0,83.517357d0,
     & 0.2511888d1,9.276449d0,82.923729d0,
     & 0.2371377d1,9.279055d0,82.592117d0,
     & 0.2238722d1,9.273203d0,82.236862d0,
     & 0.2113489d1,9.266044d0,81.857941d0,
     & 0.1995263d1,9.257412d0,81.455658d0,
     & 0.1883650d1,9.247128d0,81.030617d0,
     & 0.1778279d1,9.235045d0,80.583801d0,
     & 0.1678806d1,9.221028d0,80.116539d0,
     & 0.1584893d1,9.197192d0,79.630501d0,
     & 0.1496236d1,9.174915d0,79.151718d0,
     & 0.1412538d1,9.158957d0,78.655190d0,
     & 0.1333521d1,9.143701d0,78.143341d0,
     & 0.1258927d1,9.119092d0,77.618797d0,
     & 0.1188502d1,9.094488d0,76.978050d0,
     & 0.1122019d1,9.069882d0,76.337303d0,
     & 0.1059254d1,9.045273d0,75.696548d0,
     & 0.1000000d1,9.020666d0,75.055801d0,
     & 0.9440617d0,8.996060d0,74.415047d0,
     & 0.8912510d0,8.971451d0,73.774300d0,
     & 0.8413955d0,8.946841d0,73.133537d0,
     & 0.7943284d0,8.922239d0,72.492790d0,
     & 0.7498942d0,8.889072d0,71.754227d0,
     & 0.7079459d0,8.855903d0,71.015663d0,
     & 0.6683444d0,8.822741d0,70.277092d0,
     & 0.6309577d0,8.789569d0,69.538528d0,
     & 0.5956623d0,8.756402d0,68.799957d0,
     & 0.5623413d0,8.723238d0,68.061394d0,
     & 0.5308846d0,8.690076d0,67.322823d0,
     & 0.5011872d0,8.659819d0,66.584259d0,
     & 0.4731517d0,8.653987d0,65.812302d0,
     & 0.4466837d0,8.662815d0,64.966003d0,
     & 0.4216965d0,8.675383d0,64.033997d0,
     & 0.3981073d0,8.691288d0,63.005001d0,
     & 0.3758374d0,8.710199d0,61.867641d0,
     & 0.3548135d0,8.731828d0,60.610569d0,
     & 0.3349657d0,8.755990d0,59.222481d0,
     & 0.3162279d0,8.782597d0,57.692001d0,
     & 0.2985383d0,8.811745d0,56.007820d0,
     & 0.2818383d0,8.843682d0,54.158581d0,
     & 0.2660725d0,8.878924d0,52.132938d0,
     & 0.2511887d0,8.918329d0,49.919579d0,
     & 0.2371376d0,8.963238d0,47.507160d0,
     & 0.2238722d0,8.972139d0,44.884319d0,
     & 0.2113489d0,9.000001d0,42.152420d0,
     & 0.1995262d0,8.981428d0,39.420509d0,
     & 0.1883649d0,8.944286d0,36.688622d0,
     & 0.1778279d0,8.852855d0,33.956718d0,
     & 0.1678806d0,8.676431d0,31.224810d0,
     & 0.1584893d0,8.407143d0,28.492920d0,
     & 0.1496236d0,7.924472d0,25.651779d0,
     & 0.1412538d0,7.218974d0,23.377211d0,
     & 0.1333521d0,6.608272d0,21.564310d0,
     & 0.1258925d0,6.085515d0,20.126560d0,
     & 0.1188502d0,5.643207d0,18.992720d0,
     & 0.1122018d0,5.273376d0,18.104280d0,
     & 0.1059253d0,4.967891d0,17.413271d0,
     & 0.1000000d0,4.629994d0,16.880390d0,
     & 0.9440609d-1,4.404898d0,16.600660d0,
     & 0.8912510d-1,4.341545d0,16.353270d0,
     & 0.8413950d-1,4.290949d0,16.128790d0,
     & 0.7943282d-1,4.250250d0,15.919700d0,
     & 0.7498942d-1,4.217101d0,15.720170d0/
      data ((table(iq,ir),iq=1,3),ir=199,300)/
     & 0.7079458d-1,4.189602d0,15.525600d0,
     & 0.6683440d-1,4.166118d0,15.332490d0,
     & 0.6309573d-1,4.145311d0,15.138250d0,
     & 0.5956621d-1,4.126081d0,14.941070d0,
     & 0.5623413d-1,4.107556d0,14.739760d0,
     & 0.5308845d-1,4.088943d0,14.533670d0,
     & 0.5011873d-1,4.069606d0,14.322700d0,
     & 0.4731512d-1,4.049028d0,14.107130d0,
     & 0.4466836d-1,4.026782d0,13.887660d0,
     & 0.4216965d-1,4.002524d0,13.665300d0,
     & 0.3981072d-1,3.975116d0,13.441370d0,
     & 0.3548133d-1,3.914132d0,12.995180d0,
     & 0.3162277d-1,3.844031d0,12.563350d0,
     & 0.2818383d-1,3.765947d0,12.161090d0,
     & 0.2511886d-1,3.682950d0,11.802320d0,
     & 0.2238721d-1,3.599389d0,11.497680d0,
     & 0.1995262d-1,3.520215d0,11.252810d0,
     & 0.1778279d-1,3.450083d0,11.067720d0,
     & 0.1584893d-1,3.392399d0,10.937060d0,
     & 0.1412538d-1,3.348625d0,10.851600d0,
     & 0.1258925d-1,3.318210d0,10.800260d0,
     & 0.1122018d-1,3.299013d0,10.772190d0,
     & 0.1000000d-1,3.288110d0,10.758350d0,
     & 0.8912510d-2,3.282602d0,10.752210d0,
     & 0.7943284d-2,3.280155d0,10.749700d0,
     & 0.7079456d-2,3.279223d0,10.748640d0,
     & 0.6309573d-2,3.278931d0,10.748020d0,
     & 0.5623413d-2,3.278870d0,10.747450d0,
     & 0.5011873d-2,3.278882d0,10.746770d0,
     & 0.4466837d-2,3.278913d0,10.745920d0,
     & 0.3981071d-2,3.278958d0,10.744850d0,
     & 0.3548133d-2,3.279345d0,10.743490d0,
     & 0.3162278d-2,3.279804d0,10.741020d0,
     & 0.2818383d-2,3.280047d0,10.737780d0,
     & 0.2511887d-2,3.280304d0,10.733660d0,
     & 0.2238721d-2,3.280641d0,10.728430d0,
     & 0.1995262d-2,3.282239d0,10.721750d0,
     & 0.1778279d-2,3.283977d0,10.710610d0,
     & 0.1584893d-2,3.285009d0,10.696460d0,
     & 0.1412538d-2,3.286310d0,10.678480d0,
     & 0.1258925d-2,3.287937d0,10.655650d0,
     & 0.1122018d-2,3.289977d0,10.626710d0,
     & 0.1000000d-2,3.292506d0,10.590050d0,
     & 0.8912510d-3,3.295610d0,10.543740d0,
     & 0.7943284d-3,3.299383d0,10.485400d0,
     & 0.7079456d-3,3.303900d0,10.412190d0,
     & 0.6309572d-3,3.309220d0,10.320760d0,
     & 0.5623413d-3,3.315323d0,10.207230d0,
     & 0.5011873d-3,3.322111d0,10.067280d0,
     & 0.4466837d-3,3.329330d0,9.896221d0,
     & 0.3981071d-3,3.336508d0,9.689258d0,
     & 0.3548134d-3,3.342868d0,9.441833d0,
     & 0.3162278d-3,3.347242d0,9.150148d0,
     & 0.2818383d-3,3.347975d0,8.811836d0,
     & 0.2511887d-3,3.342849d0,8.426746d0,
     & 0.2238721d-3,3.329062d0,7.997792d0,
     & 0.1995262d-3,3.303299d0,7.531672d0,
     & 0.1778279d-3,3.261962d0,7.039310d0,
     & 0.1584893d-3,3.201615d0,6.535744d0,
     & 0.1412538d-3,3.119701d0,6.039272d0,
     & 0.1258925d-3,3.015544d0,5.569750d0,
     & 0.1122018d-3,2.891436d0,5.146163d0,
     & 0.1000000d-3,2.753449d0,4.783859d0,
     & 0.8912505d-4,2.611334d0,4.492088d0,
     & 0.7943284d-4,2.476926d0,4.272558d0,
     & 0.7079457d-4,2.361244d0,4.119564d0,
     & 0.6309577d-4,2.271389d0,4.021753d0,
     & 0.5623413d-4,2.208873d0,3.965008d0,
     & 0.5011870d-4,2.170196d0,3.935492d0,
     & 0.4466837d-4,2.149093d0,3.921912d0,
     & 0.3981071d-4,2.139039d0,3.916468d0,
     & 0.3548135d-4,2.134906d0,3.914598d0,
     & 0.3162278d-4,2.133462d0,3.914058d0,
     & 0.2818382d-4,2.133043d0,3.913929d0,
     & 0.2511887d-4,2.132938d0,3.913901d0,
     & 0.2238721d-4,2.132916d0,3.913901d0,
     & 0.1995263d-4,2.132916d0,3.913901d0,
     & 0.0000000d0,2.132916d0,3.913901d0,
     & 72*0.0d0/
      nf=276
      do i=1,300
         tgev(i)=table(1,i)
         fg(i)=table(2,i)
         fh(i)=table(3,i)
      enddo
      rdinit = 1234
      end
