from sys import argv, stdout, stderr
import ROOT
import sys
import math
import array
import numpy

def startUnc(channel):
	if channel == "gg1":
        	channelName = "boost"
		ggScale = 1.10
        else:
        	channelName = channel
		if channel == "vbf":
			ggScale = 1.30
		else:
			ggScale = 1.10
        f = open(channelName+'mutau/unc.vals','w')
	f.write(channelName + 'mutau LFVGG,SMGG126,WWGG126 Theo_PDF_gg 1.08\n')
	f.write(channelName + 'mutau LFVVBF,SMVBF126,WWVBF126 Theo_PDF_vbf 1.08\n')
	f.write('\n')
	f.write(channelName + 'mutau LFVGG, LFVVBF, SMGG126, WWGG126, SMVBF126, WWVBF126 Theo_UE 1.04\n')
	f.write('\n')
	f.write(channelName + 'mutau LFVGG,SMGG126, WWGG126 Theo_Scale_gg ' + str(ggScale)+'\n')
	f.write('\n')
	f.write(channelName + 'mutau LFVVBF,SMVBF126, WWVBF126 Theo_Scale_vbf 1.04\n')
	f.write(channelName + 'mutau signal,ww,ttbar,singlet,ztautau,zjetsother,SMVBF126,SMGG126,WWGG126,WWVBF126  lumi 1.026\n') 
	f.write(channelName + 'mutau signal,ww,ttbar,singlet,ztautau,zjetsother,SMVBF126,SMGG126,WWGG126,WWVBF126  Effi_Mu  1.02\n')
	f.write(channelName + 'mutau signal,ww,ttbar,singlet,ztautau,zjetsother,SMVBF126,SMGG126,WWGG126,WWVBF126  Effi_Tau 1.06\n') 
	f.write('\n')
	f.write('\n')
	f.write('\n')
	f.write(channelName + 'mutau ww        Norm_WW       1.15\n')
	f.write(channelName + 'mutau singlet       Norm_TOP      1.10\n')
	f.write(channelName + 'mutau ttbar        Norm_TT       1.10\n')
	f.write(channelName + 'mutau ztautau   Norm_ZTauTau  1.03\n')
	f.write(channelName + 'mutau zjetsother   Norm_ZJetsOther  1.1\n')
	f.write(channelName + 'mutau fakes     Norm_FAKES  1.30\n')
	f.write('\n')
	f.write('\n')
	f.write('\n')
	f.write(channelName + 'mutau signal,ww,singlet,ttbar,ztautau,zjetsother,SMVBF126,SMGG126 shape_JES 1\n')
	f.write(channelName + 'mutau signal,ww,singlet,ttbar,ztautau,zjetsother, SMVBF126,SMGG126 shape_TES 1\n')
        f.write(channelName + 'mutau signal,ww,singlet,ttbar,ztautau,zjetsother, SMVBF126,SMGG126 shape_MET 1\n')
	f.write('\n')
	f.write('\n')

	u = open(channelName+'mutau/unc.conf','w')
	u.write('lumi lnN\n')
	u.write('Effi_Mu lnN\n')
	u.write('Effi_Tau lnN\n')
	u.write('\n')
	u.write('Theo_UE lnN\n')
	u.write('\n')
	u.write('Theo_Scale_gg lnN\n')
	u.write('Theo_PDF_gg lnN\n')
	u.write('\n')
        u.write('Theo_Scale_vbf lnN\n')
        u.write('Theo_PDF_vbf lnN\n')
        u.write('\n')
	u.write('Norm_WW lnN\n')
	u.write('Norm_TT lnN\n')
	u.write('Norm_ZTauTau lnN\n')
	u.write('Norm_ZJetsOther lnN\n')
	u.write('Norm_FAKES lnN\n')
	u.write('Norm_SMGG lnN\n')
	u.write('Norm_SMVBF lnN\n')
        u.write('Norm_WWGG lnN\n')
        u.write('Norm_WWVBF lnN\n')
	u.write('Norm_TOP lnN\n')
	u.write('\n')
	u.write('shape_JES lnN\n')
	u.write('shape_TES lnN\n')
	u.write('shape_MET lnN\n')
	u.write('\n')
	u.write('shape_JES shape\n')
	u.write('shape_TES shape\n')
	u.write('shape_MET shape\n')
	u.write('\n')
	
	
	
def makeBinShape(histo,histoname,outfile,channel):
        maxBin = histo.GetNbinsX()
        for i in range (1,maxBin):
                if histo.GetBinContent(i) !=0:
                        #print histoname
                        #print i
                        histoBin_up = histo.Clone()
                        histoBin_down = histo.Clone()
                        histoBin_up.SetBinContent(i,(histo.GetBinContent(i)+histo.GetBinError(i)))
                        histoBin_down.SetBinContent(i,(histo.GetBinContent(i)-histo.GetBinError(i)))
			if channel == "vbf":
				channelLabel = "muhad_2Jets"
			elif channel == "gg1":
				channelLabel = "muhad_1Jet"
			elif channel == "gg0":
				channelLabel = "muhad_0Jet"
                        binName = "Stat_"+channelLabel+"_"+histoname+"_bin"+str(i)
                        histoBin_up.Write(histoname+"_"+binName+"Up")
                        histoBin_down.Write(histoname+"_"+binName+"Down")
                        if channel == "gg1":
                                channelName = "boost"
                        else:
                                channelName = channel
			print channelName+'mutau/unc.conf'
			f = open(channelName+'mutau/unc.conf','a')
			f.write(binName+" shape\n")
			f.close()
			u= open(channelName+'mutau/unc.vals','a')
			u.write(channelName+'mutau '+histoname+' '+binName+' 1\n')
			u.close()

def writeUpDownNorm(histo,histoUp,histoDown,histoname,channel,sys):
	if histo.Integral() == 0:
		return
	shiftUp = (histoUp.Integral()-histo.Integral())/histo.Integral()+1
	shiftDown = (histoDown.Integral()-histo.Integral())/histo.Integral()+1
        if channel == "gg1":
        	channelName = "boost"
        else:
        	channelName = channel
	u = open(channelName+'mutau/unc.vals','a')
	if sys == "jes":
		normName = 'shape_JES'
	elif sys == "tes":
		normName = 'shape_TES'
	elif sys == "ues":
		normName = 'shape_MET'
	else:
		print "Error: Unknown systematic"
		sys.exit()
	diffUp = histoUp.Integral()/histo.Integral()
	diffUp = abs(1-diffUp) + 1
	print diffUp
        diffDown = histoDown.Integral()/histo.Integral()
	diffDown = abs(1-diffDown) + 1
	print diffDown
	diff = max(diffUp,diffDown)
	strDiff = channelName+'mutau '+histoname+' '+normName + ' %f \n'%(diff)
	u.write(strDiff)
	

#nonedir = "jesnonesignalFeb26_01JetFix/"
#nonedir = "signal_March31_SOverSPlusB/"
#nonedir = "signal_March31_BR1/"
nonedir = "signal_March31_retest/"
jesupdir = "jesplussignalFeb26_01JetFix/"
jesdowndir = "jesminussignalFeb26_01JetFix/"
uesupdir = "uesplussignalApr26/"
uesdowndir = "uesminussignalApr26/"
tesupdir = "tesup_signal_March9/"
tesdowndir = "tesdown_signal_March9/"

channel = "vbf"

startUnc(channel)

nonefile_str = "LFV_" + channel+"_collMass_type1_fakeRate_zjetsEmbed_newSignal.root"
jesupfile_str = "LFV_" + channel+"_collMass_jes_plus_fakeRate_zjetsEmbed_newSignal_jesplus.root"
jesdownfile_str = "LFV_" + channel+"_collMass_jes_minus_fakeRate_zjetsEmbed_newSignal_jesminus.root"
uesupfile_str = "LFV_" + channel+"_collMass_type1_ues_plus_fakeRate_zjetsEmbed_newSignal.root"
uesdownfile_str = "LFV_" + channel+"_collMass_type1_ues_minus_fakeRate_zjetsEmbed_newSignal.root"
tesupfile_str = "LFV_" + channel+"_collMass_type1_fakeRate_zjetsEmbed_newSignal.root"
tesdownfile_str = "LFV_" + channel+"_collMass_type1_fakeRate_zjetsEmbed_newSignal.root"

nonefile = ROOT.TFile(nonedir+nonefile_str)
jesupfile = ROOT.TFile(jesupdir+jesupfile_str)
jesdownfile = ROOT.TFile(jesdowndir+jesdownfile_str)
uesupfile = ROOT.TFile(uesupdir+uesupfile_str)
uesdownfile = ROOT.TFile(uesdowndir+uesdownfile_str)
tesupfile = ROOT.TFile(tesupdir+tesupfile_str)
tesdownfile = ROOT.TFile(tesdowndir+tesdownfile_str)

if channel == "vbf":
	dirname = "vbfmutau"
elif channel == "gg1":
	dirname = "boostmutau"
elif channel == "gg0":
	dirname = "gg0mutau"

fakes = nonefile.Get("vbfmutau/fakes")
ztautau = nonefile.Get("vbfmutau/ztautau")
zjetsother = nonefile.Get("vbfmutau/zjetsother")
ttbar = nonefile.Get("vbfmutau/ttbar")
ww = nonefile.Get("vbfmutau/ww")
singlet = nonefile.Get("vbfmutau/singlet")
data_obs = nonefile.Get("vbfmutau/data_obs")
LFVVBF126 = nonefile.Get("vbfmutau/LFVVBF126")
LFVGG126 = nonefile.Get("vbfmutau/LFVGG126")
SMVBF126 = nonefile.Get("vbfmutau/SMVBF126")
SMGG126 = nonefile.Get("vbfmutau/SMGG126")
WWVBF126 = nonefile.Get("vbfmutau/WWVBF126")
WWGG126 = nonefile.Get("vbfmutau/WWGG126")
LFVVBF = nonefile.Get("vbfmutau/LFVVBF")
LFVGG = nonefile.Get("vbfmutau/LFVGG")
SMVBF = nonefile.Get("vbfmutau/SMVBF")
SMGG = nonefile.Get("vbfmutau/SMGG")
WWVBF = nonefile.Get("vbfmutau/WWVBF")
WWGG = nonefile.Get("vbfmutau/WWGG")

ztautau_shape_JES_Up = jesupfile.Get("vbfmutau/ztautau")
zjetsother_shape_JES_Up = jesupfile.Get("vbfmutau/zjetsother")
ttbar_shape_JES_Up = jesupfile.Get("vbfmutau/ttbar")
ww_shape_JES_Up = jesupfile.Get("vbfmutau/ww")
singlet_shape_JES_Up = jesupfile.Get("vbfmutau/singlet")
LFVVBF126_shape_JES_Up = jesupfile.Get("vbfmutau/LFVVBF126")
LFVGG126_shape_JES_Up = jesupfile.Get("vbfmutau/LFVGG126")
SMVBF126_shape_JES_Up = jesupfile.Get("vbfmutau/SMVBF126")
SMGG126_shape_JES_Up = jesupfile.Get("vbfmutau/SMGG126")
LFVVBF_shape_JES_Up = jesupfile.Get("vbfmutau/LFVVBF")
LFVGG_shape_JES_Up = jesupfile.Get("vbfmutau/LFVGG")
SMVBF_shape_JES_Up = jesupfile.Get("vbfmutau/SMVBF")
SMGG_shape_JES_Up = jesupfile.Get("vbfmutau/SMGG")

ztautau_shape_JES_Down = jesdownfile.Get("vbfmutau/ztautau")
zjetsother_shape_JES_Down = jesdownfile.Get("vbfmutau/zjetsother")
ttbar_shape_JES_Down = jesdownfile.Get("vbfmutau/ttbar")
ww_shape_JES_Down = jesdownfile.Get("vbfmutau/ww")
singlet_shape_JES_Down = jesdownfile.Get("vbfmutau/singlet")
LFVVBF126_shape_JES_Down = jesdownfile.Get("vbfmutau/LFVVBF126")
LFVGG126_shape_JES_Down = jesdownfile.Get("vbfmutau/LFVGG126")
SMVBF126_shape_JES_Down = jesdownfile.Get("vbfmutau/SMVBF126")
SMGG126_shape_JES_Down = jesdownfile.Get("vbfmutau/SMGG126")
LFVVBF_shape_JES_Down = jesdownfile.Get("vbfmutau/LFVVBF")
LFVGG_shape_JES_Down = jesdownfile.Get("vbfmutau/LFVGG")
SMVBF_shape_JES_Down = jesdownfile.Get("vbfmutau/SMVBF")
SMGG_shape_JES_Down = jesdownfile.Get("vbfmutau/SMGG")

ztautau_shape_MET_Up = uesupfile.Get("vbfmutau/ztautau")
zjetsother_shape_MET_Up = uesupfile.Get("vbfmutau/zjetsother")
ttbar_shape_MET_Up = uesupfile.Get("vbfmutau/ttbar")
ww_shape_MET_Up = uesupfile.Get("vbfmutau/ww")
singlet_shape_MET_Up = uesupfile.Get("vbfmutau/singlet")
LFVVBF126_shape_MET_Up = uesupfile.Get("vbfmutau/LFVVBF126")
LFVGG126_shape_MET_Up = uesupfile.Get("vbfmutau/LFVGG126")
SMVBF126_shape_MET_Up = uesupfile.Get("vbfmutau/SMVBF126")
SMGG126_shape_MET_Up = uesupfile.Get("vbfmutau/SMGG126")
LFVVBF_shape_MET_Up = uesupfile.Get("vbfmutau/LFVVBF")
LFVGG_shape_MET_Up = uesupfile.Get("vbfmutau/LFVGG")
SMVBF_shape_MET_Up = uesupfile.Get("vbfmutau/SMVBF")
SMGG_shape_MET_Up = uesupfile.Get("vbfmutau/SMGG")

ztautau_shape_MET_Down = uesdownfile.Get("vbfmutau/ztautau")
zjetsother_shape_MET_Down = uesdownfile.Get("vbfmutau/zjetsother")
ttbar_shape_MET_Down = uesdownfile.Get("vbfmutau/ttbar")
ww_shape_MET_Down = uesdownfile.Get("vbfmutau/ww")
singlet_shape_MET_Down = uesdownfile.Get("vbfmutau/singlet")
LFVVBF126_shape_MET_Down = uesdownfile.Get("vbfmutau/LFVVBF126")
LFVGG126_shape_MET_Down = uesdownfile.Get("vbfmutau/LFVGG126")
SMVBF126_shape_MET_Down = uesdownfile.Get("vbfmutau/SMVBF126")
SMGG126_shape_MET_Down = uesdownfile.Get("vbfmutau/SMGG126")
LFVVBF_shape_MET_Down = uesdownfile.Get("vbfmutau/LFVVBF")
LFVGG_shape_MET_Down = uesdownfile.Get("vbfmutau/LFVGG")
SMVBF_shape_MET_Down = uesdownfile.Get("vbfmutau/SMVBF")
SMGG_shape_MET_Down = uesdownfile.Get("vbfmutau/SMGG")

ztautau_shape_TES_Up = tesupfile.Get("vbfmutau/ztautau")
zjetsother_shape_TES_Up = tesupfile.Get("vbfmutau/zjetsother")
ttbar_shape_TES_Up = tesupfile.Get("vbfmutau/ttbar")
ww_shape_TES_Up = tesupfile.Get("vbfmutau/ww")
singlet_shape_TES_Up = tesupfile.Get("vbfmutau/singlet")
LFVVBF126_shape_TES_Up = tesupfile.Get("vbfmutau/LFVVBF126")
LFVGG126_shape_TES_Up = tesupfile.Get("vbfmutau/LFVGG126")
SMVBF126_shape_TES_Up = tesupfile.Get("vbfmutau/SMVBF126")
SMGG126_shape_TES_Up = tesupfile.Get("vbfmutau/SMGG126")
LFVVBF_shape_TES_Up = tesupfile.Get("vbfmutau/LFVVBF")
LFVGG_shape_TES_Up = tesupfile.Get("vbfmutau/LFVGG")
SMVBF_shape_TES_Up = tesupfile.Get("vbfmutau/SMVBF")
SMGG_shape_TES_Up = tesupfile.Get("vbfmutau/SMGG")

ztautau_shape_TES_Down = tesdownfile.Get("vbfmutau/ztautau")
zjetsother_shape_TES_Down = tesdownfile.Get("vbfmutau/zjetsother")
ttbar_shape_TES_Down = tesdownfile.Get("vbfmutau/ttbar")
ww_shape_TES_Down = tesdownfile.Get("vbfmutau/ww")
singlet_shape_TES_Down = tesdownfile.Get("vbfmutau/singlet")
LFVVBF126_shape_TES_Down = tesdownfile.Get("vbfmutau/LFVVBF126")
LFVGG126_shape_TES_Down = tesdownfile.Get("vbfmutau/LFVGG126")
SMVBF126_shape_TES_Down = tesdownfile.Get("vbfmutau/SMVBF126")
SMGG126_shape_TES_Down = tesdownfile.Get("vbfmutau/SMGG126")
LFVVBF_shape_TES_Down = tesdownfile.Get("vbfmutau/LFVVBF")
LFVGG_shape_TES_Down = tesdownfile.Get("vbfmutau/LFVGG")
SMVBF_shape_TES_Down = tesdownfile.Get("vbfmutau/SMVBF")
SMGG_shape_TES_Down = tesdownfile.Get("vbfmutau/SMGG")


##create root file with yields for datacards
outfile = ROOT.TFile(dirname+"/LFV_"+channel+"_jesshape_tesshape_uesshape_binshape.root","RECREATE")
outfile.mkdir(dirname)
outfile.cd(dirname+"/")

fakes.Write("fakes")
ztautau.Write("ztautau")
zjetsother.Write("zjetsother")
ttbar.Write("ttbar")
ww.Write("ww")
singlet.Write("singlet")
data_obs.Write("data_obs")
LFVVBF126.Write("LFVVBF126")
LFVGG126.Write("LFVGG126")
SMVBF126.Write("SMVBF126")
SMGG126.Write("SMGG126")
WWVBF126.Write("WWVBF126")
WWGG126.Write("WWGG126")
LFVVBF.Write("LFVVBF")
LFVGG.Write("LFVGG")
SMVBF.Write("SMVBF")
SMGG.Write("SMGG")
WWVBF.Write("WWVBF")
WWGG.Write("WWGG")

makeBinShape(fakes,"fakes",outfile,channel)
makeBinShape(ztautau,"ztautau",outfile,channel)
makeBinShape(zjetsother,"zjetsother",outfile,channel)
makeBinShape(ttbar,"ttbar",outfile,channel)
makeBinShape(ww,"ww",outfile,channel)
makeBinShape(singlet,"singlet",outfile,channel)
makeBinShape(LFVVBF,"LFVVBF",outfile,channel)
makeBinShape(LFVGG,"LFVGG",outfile,channel)
makeBinShape(SMVBF126,"SMVBF126",outfile,channel)
makeBinShape(SMGG126,"SMGG126",outfile,channel)
makeBinShape(WWVBF126,"WWVBF126",outfile,channel)
makeBinShape(WWGG126,"WWGG126",outfile,channel)

writeUpDownNorm(ztautau,ztautau_shape_JES_Up,ztautau_shape_JES_Down,"ztautau",channel,"jes")
writeUpDownNorm(zjetsother,zjetsother_shape_JES_Up,zjetsother_shape_JES_Down,"zjetsother",channel,"jes")
writeUpDownNorm(ttbar,ttbar_shape_JES_Up,ttbar_shape_JES_Down,"ttbar",channel,"jes")
writeUpDownNorm(ww,ww_shape_JES_Up,ww_shape_JES_Down,"ww",channel,"jes")
writeUpDownNorm(singlet,singlet_shape_JES_Up,singlet_shape_JES_Down,"singlet",channel,"jes")
writeUpDownNorm(SMVBF126,SMVBF126_shape_JES_Up,SMVBF126_shape_JES_Down,"SMVBF126",channel,"jes")
writeUpDownNorm(SMGG126,SMGG126_shape_JES_Up,SMGG126_shape_JES_Down,"SMGG126",channel,"jes")
writeUpDownNorm(LFVVBF,LFVVBF_shape_JES_Up,LFVVBF_shape_JES_Down,"LVFVBF",channel,"jes")
writeUpDownNorm(LFVGG,LFVGG_shape_JES_Up,LFVGG_shape_JES_Down,"LFVGG",channel,"jes")

writeUpDownNorm(ztautau,ztautau_shape_JES_Up,ztautau_shape_JES_Down,"ztautau",channel,"ues")
writeUpDownNorm(zjetsother,zjetsother_shape_JES_Up,zjetsother_shape_JES_Down,"zjetsother",channel,"ues")
writeUpDownNorm(ttbar,ttbar_shape_JES_Up,ttbar_shape_JES_Down,"ttbar",channel,"ues")
writeUpDownNorm(ww,ww_shape_JES_Up,ww_shape_JES_Down,"ww",channel,"ues")
writeUpDownNorm(singlet,singlet_shape_JES_Up,singlet_shape_JES_Down,"singlet",channel,"ues")
writeUpDownNorm(SMVBF126,SMVBF126_shape_JES_Up,SMVBF126_shape_JES_Down,"SMVBF126",channel,"ues")
writeUpDownNorm(SMGG126,SMGG126_shape_JES_Up,SMGG126_shape_JES_Down,"SMGG126",channel,"ues")
writeUpDownNorm(LFVVBF,LFVVBF_shape_JES_Up,LFVVBF_shape_JES_Down,"LVFVBF",channel,"ues")
writeUpDownNorm(LFVGG,LFVGG_shape_JES_Up,LFVGG_shape_JES_Down,"LFVGG",channel,"ues")

writeUpDownNorm(ztautau,ztautau_shape_TES_Up,ztautau_shape_TES_Down,"ztautau",channel,"tes")
writeUpDownNorm(zjetsother,zjetsother_shape_TES_Up,zjetsother_shape_TES_Down,"zjetsother",channel,"tes")
writeUpDownNorm(ttbar,ttbar_shape_TES_Up,ttbar_shape_TES_Down,"ttbar",channel,"tes")
writeUpDownNorm(ww,ww_shape_TES_Up,ww_shape_TES_Down,"ww",channel,"tes")
writeUpDownNorm(singlet,singlet_shape_TES_Up,singlet_shape_TES_Down,"singlet",channel,"tes")
writeUpDownNorm(SMVBF126,SMVBF126_shape_TES_Up,SMVBF126_shape_TES_Down,"SMVBF126",channel,"tes")
writeUpDownNorm(SMGG126,SMGG126_shape_TES_Up,SMGG126_shape_TES_Down,"SMGG126",channel,"tes")
writeUpDownNorm(LFVVBF,LFVVBF_shape_TES_Up,LFVVBF_shape_TES_Down,"LVFVBF",channel,"tes")
writeUpDownNorm(LFVGG,LFVGG_shape_TES_Up,LFVGG_shape_TES_Down,"LFVGG",channel,"tes")
#scale up down shapes to original 


ztautau_shape_JES_Up.Scale(ztautau.Integral()/ztautau_shape_JES_Up.Integral())
zjetsother_shape_JES_Up.Scale(zjetsother.Integral()/zjetsother_shape_JES_Up.Integral())
ttbar_shape_JES_Up.Scale(ttbar.Integral()/ttbar_shape_JES_Up.Integral())
ww_shape_JES_Up.Scale(ww.Integral()/ww_shape_JES_Up.Integral())
singlet_shape_JES_Up.Scale(singlet.Integral()/singlet_shape_JES_Up.Integral())
LFVVBF126_shape_JES_Up.Scale(LFVVBF126.Integral()/LFVVBF126_shape_JES_Up.Integral())
LFVGG126_shape_JES_Up.Scale(LFVGG126.Integral()/LFVGG126_shape_JES_Up.Integral())
SMVBF126_shape_JES_Up.Scale(SMVBF126.Integral()/SMVBF126_shape_JES_Up.Integral())
SMGG126_shape_JES_Up.Scale(SMGG126.Integral()/SMGG126_shape_JES_Up.Integral())
LFVVBF_shape_JES_Up.Scale(LFVVBF.Integral()/LFVVBF_shape_JES_Up.Integral())
LFVGG_shape_JES_Up.Scale(LFVGG.Integral()/LFVGG_shape_JES_Up.Integral())
SMVBF_shape_JES_Up.Scale(SMVBF.Integral()/SMVBF_shape_JES_Up.Integral())
SMGG_shape_JES_Up.Scale(SMGG.Integral()/SMGG_shape_JES_Up.Integral())

ztautau_shape_JES_Down.Scale(ztautau.Integral()/ztautau_shape_JES_Down.Integral())
zjetsother_shape_JES_Down.Scale(zjetsother.Integral()/zjetsother_shape_JES_Down.Integral())
ttbar_shape_JES_Down.Scale(ttbar.Integral()/ttbar_shape_JES_Down.Integral())
ww_shape_JES_Down.Scale(ww.Integral()/ww_shape_JES_Down.Integral())
singlet_shape_JES_Down.Scale(singlet.Integral()/singlet_shape_JES_Down.Integral())
LFVVBF126_shape_JES_Down.Scale(LFVVBF126.Integral()/LFVVBF126_shape_JES_Down.Integral())
LFVGG126_shape_JES_Down.Scale(LFVGG126.Integral()/LFVGG126_shape_JES_Down.Integral())
SMVBF126_shape_JES_Down.Scale(SMVBF126.Integral()/SMVBF126_shape_JES_Down.Integral())
SMGG126_shape_JES_Down.Scale(SMGG126.Integral()/SMGG126_shape_JES_Down.Integral())
LFVVBF_shape_JES_Down.Scale(LFVVBF.Integral()/LFVVBF_shape_JES_Down.Integral())
LFVGG_shape_JES_Down.Scale(LFVGG.Integral()/LFVGG_shape_JES_Down.Integral())
SMVBF_shape_JES_Down.Scale(SMVBF.Integral()/SMVBF_shape_JES_Down.Integral())
SMGG_shape_JES_Down.Scale(SMGG.Integral()/SMGG_shape_JES_Down.Integral())

ztautau_shape_MET_Up.Scale(ztautau.Integral()/ztautau_shape_MET_Up.Integral())
zjetsother_shape_MET_Up.Scale(zjetsother.Integral()/zjetsother_shape_MET_Up.Integral())
ttbar_shape_MET_Up.Scale(ttbar.Integral()/ttbar_shape_MET_Up.Integral())
ww_shape_MET_Up.Scale(ww.Integral()/ww_shape_MET_Up.Integral())
singlet_shape_MET_Up.Scale(singlet.Integral()/singlet_shape_MET_Up.Integral())
LFVVBF126_shape_MET_Up.Scale(LFVVBF126.Integral()/LFVVBF126_shape_MET_Up.Integral())
LFVGG126_shape_MET_Up.Scale(LFVGG126.Integral()/LFVGG126_shape_MET_Up.Integral())
SMVBF126_shape_MET_Up.Scale(SMVBF126.Integral()/SMVBF126_shape_MET_Up.Integral())
SMGG126_shape_MET_Up.Scale(SMGG126.Integral()/SMGG126_shape_MET_Up.Integral())
LFVVBF_shape_MET_Up.Scale(LFVVBF.Integral()/LFVVBF_shape_MET_Up.Integral())
LFVGG_shape_MET_Up.Scale(LFVGG.Integral()/LFVGG_shape_MET_Up.Integral())
SMVBF_shape_MET_Up.Scale(SMVBF.Integral()/SMVBF_shape_MET_Up.Integral())
SMGG_shape_MET_Up.Scale(SMGG.Integral()/SMGG_shape_MET_Up.Integral())

ztautau_shape_MET_Down.Scale(ztautau.Integral()/ztautau_shape_MET_Down.Integral())
zjetsother_shape_MET_Down.Scale(zjetsother.Integral()/zjetsother_shape_MET_Down.Integral())
ttbar_shape_MET_Down.Scale(ttbar.Integral()/ttbar_shape_MET_Down.Integral())
ww_shape_MET_Down.Scale(ww.Integral()/ww_shape_MET_Down.Integral())
singlet_shape_MET_Down.Scale(singlet.Integral()/singlet_shape_MET_Down.Integral())
LFVVBF126_shape_MET_Down.Scale(LFVVBF126.Integral()/LFVVBF126_shape_MET_Down.Integral())
LFVGG126_shape_MET_Down.Scale(LFVGG126.Integral()/LFVGG126_shape_MET_Down.Integral())
SMVBF126_shape_MET_Down.Scale(SMVBF126.Integral()/SMVBF126_shape_MET_Down.Integral())
SMGG126_shape_MET_Down.Scale(SMGG126.Integral()/SMGG126_shape_MET_Down.Integral())
LFVVBF_shape_MET_Down.Scale(LFVVBF.Integral()/LFVVBF_shape_MET_Down.Integral())
LFVGG_shape_MET_Down.Scale(LFVGG.Integral()/LFVGG_shape_MET_Down.Integral())
SMVBF_shape_MET_Down.Scale(SMVBF.Integral()/SMVBF_shape_MET_Down.Integral())
SMGG_shape_MET_Down.Scale(SMGG.Integral()/SMGG_shape_MET_Down.Integral())

ztautau_shape_TES_Up.Scale(ztautau.Integral()/ztautau_shape_TES_Up.Integral())
zjetsother_shape_TES_Up.Scale(zjetsother.Integral()/zjetsother_shape_TES_Up.Integral())
ttbar_shape_TES_Up.Scale(ttbar.Integral()/ttbar_shape_TES_Up.Integral())
ww_shape_TES_Up.Scale(ww.Integral()/ww_shape_TES_Up.Integral())
singlet_shape_TES_Up.Scale(singlet.Integral()/singlet_shape_TES_Up.Integral())
LFVVBF126_shape_TES_Up.Scale(LFVVBF126.Integral()/LFVVBF126_shape_TES_Up.Integral())
LFVGG126_shape_TES_Up.Scale(LFVGG126.Integral()/LFVGG126_shape_TES_Up.Integral())
SMVBF126_shape_TES_Up.Scale(SMVBF126.Integral()/SMVBF126_shape_TES_Up.Integral())
SMGG126_shape_TES_Up.Scale(SMGG126.Integral()/SMGG126_shape_TES_Up.Integral())
LFVVBF_shape_TES_Up.Scale(LFVVBF.Integral()/LFVVBF_shape_TES_Up.Integral())
LFVGG_shape_TES_Up.Scale(LFVGG.Integral()/LFVGG_shape_TES_Up.Integral())
SMVBF_shape_TES_Up.Scale(SMVBF.Integral()/SMVBF_shape_TES_Up.Integral())
SMGG_shape_TES_Up.Scale(SMGG.Integral()/SMGG_shape_TES_Up.Integral())

ztautau_shape_TES_Down.Scale(ztautau.Integral()/ztautau_shape_TES_Down.Integral())
zjetsother_shape_TES_Down.Scale(zjetsother.Integral()/zjetsother_shape_TES_Down.Integral())
ttbar_shape_TES_Down.Scale(ttbar.Integral()/ttbar_shape_TES_Down.Integral())
ww_shape_TES_Down.Scale(ww.Integral()/ww_shape_TES_Down.Integral())
singlet_shape_TES_Down.Scale(singlet.Integral()/singlet_shape_TES_Down.Integral())
LFVVBF126_shape_TES_Down.Scale(LFVVBF126.Integral()/LFVVBF126_shape_TES_Down.Integral())
LFVGG126_shape_TES_Down.Scale(LFVGG126.Integral()/LFVGG126_shape_TES_Down.Integral())
SMVBF126_shape_TES_Down.Scale(SMVBF126.Integral()/SMVBF126_shape_TES_Down.Integral())
SMGG126_shape_TES_Down.Scale(SMGG126.Integral()/SMGG126_shape_TES_Down.Integral())
LFVVBF_shape_TES_Down.Scale(LFVVBF.Integral()/LFVVBF_shape_TES_Down.Integral())
LFVGG_shape_TES_Down.Scale(LFVGG.Integral()/LFVGG_shape_TES_Down.Integral())
SMVBF_shape_TES_Down.Scale(SMVBF.Integral()/SMVBF_shape_TES_Down.Integral())
SMGG_shape_TES_Down.Scale(SMGG.Integral()/SMGG_shape_TES_Down.Integral())

ztautau_shape_JES_Up.Write("ztautau_shape_JESUp")
zjetsother_shape_JES_Up.Write("zjetsother_shape_JESUp")
ttbar_shape_JES_Up.Write("ttbar_shape_JESUp")
ww_shape_JES_Up.Write("ww_shape_JESUp")
singlet_shape_JES_Up.Write("singlet_shape_JESUp")
LFVVBF126_shape_JES_Up.Write("LFVVBF126_shape_JESUp")
LFVGG126_shape_JES_Up.Write("LFVGG126_shape_JESUp")
SMVBF126_shape_JES_Up.Write("SMVBF126_shape_JESUp")
SMGG126_shape_JES_Up.Write("SMGG126_shape_JESUp")
LFVVBF_shape_JES_Up.Write("LFVVBF_shape_JESUp")
LFVGG_shape_JES_Up.Write("LFVGG_shape_JESUp")
SMVBF_shape_JES_Up.Write("SMVBF_shape_JESUp")
SMGG_shape_JES_Up.Write("SMGG_shape_JESUp")

ztautau_shape_JES_Down.Write("ztautau_shape_JESDown")
zjetsother_shape_JES_Down.Write("zjetsother_shape_JESDown")
ttbar_shape_JES_Down.Write("ttbar_shape_JESDown")
ww_shape_JES_Down.Write("ww_shape_JESDown")
singlet_shape_JES_Down.Write("singlet_shape_JESDown")
LFVVBF126_shape_JES_Down.Write("LFVVBF126_shape_JESDown")
LFVGG126_shape_JES_Down.Write("LFVGG126_shape_JESDown")
SMVBF126_shape_JES_Down.Write("SMVBF126_shape_JESDown")
SMGG126_shape_JES_Down.Write("SMGG126_shape_JESDown")
LFVVBF_shape_JES_Down.Write("LFVVBF_shape_JESDown")
LFVGG_shape_JES_Down.Write("LFVGG_shape_JESDown")
SMVBF_shape_JES_Down.Write("SMVBF_shape_JESDown")
SMGG_shape_JES_Down.Write("SMGG_shape_JESDown")

ztautau_shape_MET_Up.Write("ztautau_shape_METUp")
zjetsother_shape_MET_Up.Write("zjetsother_shape_METUp")
ttbar_shape_MET_Up.Write("ttbar_shape_METUp")
ww_shape_MET_Up.Write("ww_shape_METUp")
singlet_shape_MET_Up.Write("singlet_shape_METUp")
LFVVBF126_shape_MET_Up.Write("LFVVBF126_shape_METUp")
LFVGG126_shape_MET_Up.Write("LFVGG126_shape_METUp")
SMVBF126_shape_MET_Up.Write("SMVBF126_shape_METUp")
SMGG126_shape_MET_Up.Write("SMGG126_shape_METUp")
LFVVBF_shape_MET_Up.Write("LFVVBF_shape_METUp")
LFVGG_shape_MET_Up.Write("LFVGG_shape_METUp")
SMVBF_shape_MET_Up.Write("SMVBF_shape_METUp")
SMGG_shape_MET_Up.Write("SMGG_shape_METUp")

ztautau_shape_MET_Down.Write("ztautau_shape_METDown")
zjetsother_shape_MET_Down.Write("zjetsother_shape_METDown")
ttbar_shape_MET_Down.Write("ttbar_shape_METDown")
ww_shape_MET_Down.Write("ww_shape_METDown")
singlet_shape_MET_Down.Write("singlet_shape_METDown")
LFVVBF126_shape_MET_Down.Write("LFVVBF126_shape_METDown")
LFVGG126_shape_MET_Down.Write("LFVGG126_shape_METDown")
SMVBF126_shape_MET_Down.Write("SMVBF126_shape_METDown")
SMGG126_shape_MET_Down.Write("SMGG126_shape_METDown")
LFVVBF_shape_MET_Down.Write("LFVVBF_shape_METDown")
LFVGG_shape_MET_Down.Write("LFVGG_shape_METDown")
SMVBF_shape_MET_Down.Write("SMVBF_shape_METDown")
SMGG_shape_MET_Down.Write("SMGG_shape_METDown")

ztautau_shape_TES_Up.Write("ztautau_shape_TESUp")
zjetsother_shape_TES_Up.Write("zjetsother_shape_TESUp")
ttbar_shape_TES_Up.Write("ttbar_shape_TESUp")
ww_shape_TES_Up.Write("ww_shape_TESUp")
singlet_shape_TES_Up.Write("singlet_shape_TESUp")
LFVVBF126_shape_TES_Up.Write("LFVVBF126_shape_TESUp")
LFVGG126_shape_TES_Up.Write("LFVGG126_shape_TESUp")
SMVBF126_shape_TES_Up.Write("SMVBF126_shape_TESUp")
SMGG126_shape_TES_Up.Write("SMGG126_shape_TESUp")
LFVVBF_shape_TES_Up.Write("LFVVBF_shape_TESUp")
LFVGG_shape_TES_Up.Write("LFVGG_shape_TESUp")
SMVBF_shape_TES_Up.Write("SMVBF_shape_TESUp")
SMGG_shape_TES_Up.Write("SMGG_shape_TESUp")

ztautau_shape_TES_Down.Write("ztautau_shape_TESDown")
zjetsother_shape_TES_Down.Write("zjetsother_shape_TESDown")
ttbar_shape_TES_Down.Write("ttbar_shape_TESDown")
ww_shape_TES_Down.Write("ww_shape_TESDown")
singlet_shape_TES_Down.Write("singlet_shape_TESDown")
LFVVBF126_shape_TES_Down.Write("LFVVBF126_shape_TESDown")
LFVGG126_shape_TES_Down.Write("LFVGG126_shape_TESDown")
SMVBF126_shape_TES_Down.Write("SMVBF126_shape_TESDown")
SMGG126_shape_TES_Down.Write("SMGG126_shape_TESDown")
LFVVBF_shape_TES_Down.Write("LFVVBF_shape_TESDown")
LFVGG_shape_TES_Down.Write("LFVGG_shape_TESDown")
SMVBF_shape_TES_Down.Write("SMVBF_shape_TESDown")
SMGG_shape_TES_Down.Write("SMGG_shape_TESDown")


outfile.Write()

