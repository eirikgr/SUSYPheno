# LaTeX2HTML site specific configuration file
# generated by config.pl

# You may edit this file to get around deficiencies of the configuration
# procedure, but you have to be sure of what you are doing!
# If you think there are bugs in the configuration procedure, please report
# them. See the BUGS file on how to do it. Your help is appreciated!

package cfgcache;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(%cfg);

$cfg{'ANYTOPNM'} = q'/usr/local/bin/anytopnm';
$cfg{'BINDIR'} = q'/usr/local/bin';
$cfg{'BMPTOPPM'} = q'/usr/local/bin/bmptoppm';
$cfg{'CRAYOLAFILE'} = q'/usr/local/share/lib/latex2html/styles/crayola.txt';
$cfg{'DVIPS'} = q'/usr/local/teTeX/bin/powerpc-apple-darwin-current/dvips';
$cfg{'DVIPSOPT'} = q' -Ppdf  -E';
$cfg{'GIFTOPNM'} = q'/usr/local/bin/giftopnm';
$cfg{'GS'} = q'/usr/local/bin/gs';
$cfg{'GSALIASDEVICE'} = q'ppmraw';
$cfg{'GSDEVICE'} = q'ppmraw';
$cfg{'GSLANDSCAPE'} = q'';
$cfg{'GS_LIB'} = q'';
$cfg{'HASHBANG'} = q'1';
$cfg{'HTML_VALIDATOR'} = q'';
$cfg{'ICONPATH'} = q'file:/usr/local/share/lib/latex2html/icons';
$cfg{'ICONSERVER'} = q'';
$cfg{'ICONSTORAGE'} = q'';
$cfg{'IMAGE_TYPES'} = q'png gif';
$cfg{'INITEX'} = q'';
$cfg{'JPEGTOPNM'} = q'/usr/local/bin/jpegtopnm';
$cfg{'KPSEWHICH'} = q'/usr/local/teTeX/bin/powerpc-apple-darwin-current/kpsewhich';
$cfg{'LATEX'} = q'/usr/local/teTeX/bin/powerpc-apple-darwin-current/latex';
$cfg{'LATEX2HTMLDIR'} = q'/usr/local/share/lib/latex2html';
$cfg{'LATEX2HTMLPLATDIR'} = q'/usr/local/lib/latex2html';
$cfg{'LIBDIR'} = q'/usr/local/lib/latex2html';
$cfg{'METADPI'} = q'0';
$cfg{'METAMODE'} = q'';
$cfg{'MKTEXLSR'} = q'/usr/local/teTeX/bin/powerpc-apple-darwin-current/mktexlsr';
$cfg{'NULLFILE'} = q'/dev/null';
$cfg{'PBMMAKE'} = q'/usr/local/bin/pbmmake';
$cfg{'PCXTOPPM'} = q'/usr/local/bin/pcxtoppm';
$cfg{'PERL'} = q'/usr/bin/perl';
$cfg{'PERLFOOTER'} = q'';
$cfg{'PERLHEADER'} = <<'EOQ';
#! /usr/bin/perl -w
EOQ
$cfg{'PERLSCRIPTDIR'} = q'/usr/bin';
$cfg{'PICTTOPPM'} = q'/usr/local/bin/picttoppm';
$cfg{'PK_GENERATION'} = q'0';
$cfg{'PNGTOPNM'} = q'/usr/local/bin/pngtopnm';
$cfg{'PNMBLACK'} = q'';
$cfg{'PNMCAT'} = q'/usr/local/bin/pnmcat';
$cfg{'PNMCROP'} = q'/usr/local/bin/pnmcrop -verbose ';
$cfg{'PNMCROPOPT'} = q' -sides ';
$cfg{'PNMCUT'} = q'/usr/local/bin/pnmcut';
$cfg{'PNMFILE'} = q'/usr/local/bin/pnmfile';
$cfg{'PNMFLIP'} = q'/usr/local/bin/pnmflip';
$cfg{'PNMPAD'} = q'/usr/local/bin/pnmpad';
$cfg{'PNMROTATE'} = q'/usr/local/bin/pnmrotate';
$cfg{'PNMSCALE'} = q'/usr/local/bin/pnmscale';
$cfg{'PNMTOPNG'} = q'/usr/local/bin/pnmtopng';
$cfg{'PPMQUANT'} = q'/usr/local/bin/ppmquant';
$cfg{'PPMTOGIF'} = q'/usr/local/bin/ppmtogif';
$cfg{'PPMTOJPEG'} = q'/usr/local/bin/ppmtojpeg';
$cfg{'PREFIX'} = q'/usr/local';
$cfg{'RGBCOLORFILE'} = q'/usr/local/share/lib/latex2html/styles/rgb.txt';
$cfg{'SGITOPNM'} = q'/usr/local/bin/sgitopnm';
$cfg{'SHLIBDIR'} = q'/usr/local/share/lib/latex2html';
$cfg{'TEX'} = q'/usr/local/teTeX/bin/powerpc-apple-darwin-current/tex';
$cfg{'TEXPATH'} = q'/usr/local/teTeX/share/texmf.local/tex/latex/html';
$cfg{'TIFFTOPNM'} = q'/usr/local/bin/tifftopnm';
$cfg{'TMPSPACE'} = q'/tmp';
$cfg{'WEB2C'} = q'1';
$cfg{'XBMTOPBM'} = q'/usr/local/bin/xbmtopbm';
$cfg{'XWDTOPNM'} = q'/usr/local/bin/xwdtopnm';
$cfg{'dd'} = q'/';
$cfg{'distver'} = q'2002-2-1';
$cfg{'exec_extension'} = q'';
$cfg{'gif_interlace'} = q'netpbm';
$cfg{'gif_trans'} = q'netpbm';
$cfg{'have_dvipsmode'} = q'';
$cfg{'have_geometry'} = q'1';
$cfg{'have_images'} = q'1';
$cfg{'have_pstoimg'} = q'1';
$cfg{'perl_starter'} = q'';
$cfg{'pipes'} = q'1';
$cfg{'plat'} = q'unix';
$cfg{'scriptdir'} = q'/usr/local/bin';
$cfg{'scriptext'} = q'';
$cfg{'srcdir'} = q'/Volumes/home/edsjo/DarkSUSY/library/trunk/contrib/latex2html-2002-2-1';
$cfg{'texlive'} = q'0';
$cfg{'wrapper'} = q'0';

1; # must be last line