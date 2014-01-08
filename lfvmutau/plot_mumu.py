from sys import argv, stdout, stderr
import ROOT
import sys
import array

ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


#Directory Structure:
#all mmt root files within mmtdir
#preselection files in preseldir
#signal files in signaldir
#files with isolation applied in isodir
#files without isolation in noisodir

#mmtdir = "MuMuTau_Sept26/mmt/"
#mmtdir = "MuMuTau_Sept27_mPt20/mmt/"
#mmtdir = "MuMuTau_Sept27_mPt20_mu17mu8/mmt/"
mmtdir = "MuMuTau_Sept30_loosevbf/mmt/"
#mmtdir = "MuMuTau_NoVBF_Old/mmt/"
#mmtdir = "MuMuTau_Oct29_twojetspresel/mmt/"
preseldir = mmtdir + "preselection/"
signaldir = mmtdir + "signal/"
preselisodir = preseldir+"isotrue/AnalyzeMuMuTauTight/"
preselnoisodir = preseldir +"isofalse/AnalyzeMuMuTauTight/"
signalisodir = signaldir+"isotrue/AnalyzeMuMuTauTight/"
signalnoisodir = signaldir +"isofalse/AnalyzeMuMuTauTight/"

def make_fakerate(isodir, noisodir, selectiondir):
	data_file_str_iso = isodir+"datammt_2012.root"
	z_file_str_iso = isodir+"Zjetsmmtvbf.root"
	data_file_str_noiso = noisodir+"datammt_2012.root"
	z_file_str_noiso = noisodir+"Zjetsmmtvbf.root"

	data_file_idiso = ROOT.TFile(data_file_str_iso)
	data_file_id = ROOT.TFile(data_file_str_noiso)
	z_file_idiso = ROOT.TFile(z_file_str_iso)
	z_file_id = ROOT.TFile(z_file_str_noiso)
	var = "tJetPt"
	doLog = True
	binwidth = 1
	channel = "vbf"

	data_idiso = data_file_idiso.Get(channel+"/"+var).Clone()
	data_id = data_file_id.Get(channel+"/"+var).Clone()
	z_idiso = z_file_idiso.Get(channel+"/"+var).Clone()
	z_id = z_file_id.Get(channel+"/"+var).Clone() 

	data1_lumifile = 'lumicalc_new/data_SingleMu_Run2012A_22Jan2013_v1_2.lumicalc.sum'
	data2_lumifile = 'lumicalc_new/data_SingleMu_Run2012B_22Jan2013_v1_2.lumicalc.sum'
	data3_lumifile = 'lumicalc_new/data_SingleMu_Run2012C_22Jan2013_v1_2.lumicalc.sum'
	data4_lumifile = 'lumicalc_new/data_SingleMu_Run2012D_22Jan2013_v1.lumicalc.sum'
	zjets_lumifile = 'lumicalc_new/Zjets_M50.lumicalc.sum'

	f = open(zjets_lumifile).read().splitlines()
	zjets_efflumi = float(f[0])

	f = open(data1_lumifile).read().splitlines()
	data1_lumi = float(f[0])

	f = open(data2_lumifile).read().splitlines()
	data2_lumi = float(f[0])

	f = open(data3_lumifile).read().splitlines()
	data3_lumi = float(f[0])

	f = open(data4_lumifile).read().splitlines()
	data4_lumi = float(f[0])

	lumi = data1_lumi+data2_lumi+data3_lumi+data4_lumi

	zjets_norm = lumi/zjets_efflumi

	z_idiso.Scale(zjets_norm)
	z_id.Scale(zjets_norm)

        canvasid = ROOT.TCanvas("canvasid","canvasid",800,800)
        canvasidiso = ROOT.TCanvas("canvasidiso","canvasidiso",800,800)
        canvasftau = ROOT.TCanvas("canvasftau","canvasftau",800,800)
        if doLog == True:
                canvasid.SetLogy()
                canvasidiso.SetLogy()
                canvasftau.SetLogy()

	z_idiso.SetFillColor(ROOT.kGreen+3)
 	z_id.SetFillColor(ROOT.kGreen+3)

	z_idiso.SetMarkerSize(0)
	z_id.SetMarkerSize(0)
        #z_idiso.Rebin(binwidth)
        #z_id.Rebin(binwidth)
        #data_idiso.Rebin(binwidth)
        #data_id.Rebin(binwidth)
	#z_ftau = z_idiso.Clone()
	#z_ftau.Divide(z_id)
	#data_ftau = data_idiso.Clone()
	#data_ftau.Divide(data_id)
	lowBinBound = 60
	medBinBound = 80
	i = 0
	bincount = -1
	binning = array.array('d',[])
	#rebin to prevent low statistics problems at high tjetpt
	while (i <= 480):
		if i < 70:
			binning.append(i)
			i = i+2
		elif i< 100:
			binning.append(i)
			i = i+5
		elif i<120:
			binning.append(i)
			i = i+10
		else:
			binning.append(i)
			i = i+20
		bincount = bincount+1
	print bincount
	print len(binning)
	print binning
	z_idisoRebin = z_idiso.Rebin(bincount,"z_idisoRebin",binning)
        z_idRebin = z_id.Rebin(bincount,"z_idRebin",binning)
	data_idisoRebin = data_idiso.Rebin(bincount,"data_idisoRebin",binning)
	data_idRebin = data_id.Rebin(bincount,"data_idRein",binning)
	z_ftau = z_idisoRebin.Clone()
	z_ftau.Divide(z_idRebin)
	data_ftau = data_idisoRebin.Clone()
	data_ftau.Divide(data_idRebin)
	legendid = ROOT.TLegend(0.45,0.6,0.8,0.8,'','brNDC')
	legendid.AddEntry(z_idRebin,'Z Mu Mu')
	legendid.AddEntry(data_idRebin,'Data')
	legendid.SetFillColor(0)
	legendid.SetBorderSize(0)
        legendidiso = ROOT.TLegend(0.45,0.6,0.8,0.8,'','brNDC')
        legendidiso.AddEntry(z_idisoRebin,'Z Mu Mu')
        legendidiso.AddEntry(data_idisoRebin,'Data')
        legendidiso.SetFillColor(0)
        legendidiso.SetBorderSize(0)
        legendftau = ROOT.TLegend(0.45,0.6,0.8,0.8,'','brNDC')
        legendftau.AddEntry(z_ftau,'Z Mu Mu')
        legendftau.AddEntry(data_ftau,'Data')
        legendftau.SetFillColor(0)
        legendftau.SetBorderSize(0)
	canvasid.cd()
	z_idRebin.Draw('hist')
	maxZ = z_idRebin.GetMaximum()
	maxData=data_idRebin.GetMaximum()
	if maxData>maxZ:
		maxHist = maxData
	else:
		maxHist = maxZ
	z_idRebin.SetMaximum(maxHist*1.05)
	z_idRebin.GetXaxis().SetTitle("Tau Jet P_{T} (GeV)")
	z_idRebin.GetYaxis().SetTitleOffset(1.2)
	lumifb = '%.2f'%(lumi/1000)
	plot_str_id = "ID"
        title_str_id = "\sqrt{s} = 8 TeV   L = " + str(lumifb)+" fb^{-1}"
        z_idRebin.SetTitle("")
        titleText = ROOT.TPaveText(0.2,0.91,0.7,0.99,"brNDC")
        titleText.AddText(title_str_id)
        titleText.SetFillStyle(0)
        titleText.Draw('sames')
	data_idRebin.Draw('sames')
	legendid.Draw('sames')
	savestr = "ID.png"
	if doLog == True:
		savestr = "Log"+savestr
	canvasid.SaveAs(noisodir+savestr)
        canvasidiso.cd()
        z_idisoRebin.Draw('hist')
        maxZ = z_idisoRebin.GetMaximum()
        maxData=data_idisoRebin.GetMaximum()
        if maxData>maxZ:
                maxHist = maxData
        else:
                maxHist = maxZ
        z_idisoRebin.SetMaximum(maxHist*1.05)
        z_idisoRebin.GetXaxis().SetTitle("Tau Jet P_{T} (GeV)")
        z_idisoRebin.GetYaxis().SetTitleOffset(1.2)
        lumifb = '%.2f'%(lumi/1000)
        plot_str_idiso = "ID+ISO"
        z_idisoRebin.SetTitle("")
        titleText.Draw('sames')
        data_idisoRebin.Draw('sames')
        legendidiso.Draw('sames')
        savestr = "IDISO.png"
        if doLog == True:
                savestr = "Log"+savestr
        canvasidiso.SaveAs(isodir+savestr)
        canvasftau.cd()
        z_ftau.Draw('hist')
        maxZ = z_ftau.GetMaximum()
        maxData=data_ftau.GetMaximum()
        if maxData>maxZ:
                maxHist = maxData
        else:
                maxHist = maxZ
        z_ftau.SetMaximum(maxHist*1.05)
        z_ftau.GetXaxis().SetTitle("Tau Jet P_{T} (GeV)")
	z_ftau.GetXaxis().SetRangeUser(0,200)
        z_ftau.GetYaxis().SetTitleOffset(1.1)
	z_ftau.GetXaxis().SetLabelSize(0.04)
	canvasftau.SetLeftMargin(0.15)
	#z_ftau.GetYaxis().SetTitleSize(0.05)
	ftauYTitle = "Tau Fake Rate"
	z_ftau.GetYaxis().SetTitle(ftauYTitle)
	#for i in range(1,200):
		#/print (z_ftau.GetBinContent(i)-z_ftau.GetBinError(i))
	
        lumifb = '%.2f'%(lumi/1000)
        plot_str_ftau = "Tau Fake Rate"
        title_str_ftau = plot_str_ftau+"\sqrt{8} TeV Collisions  L = " + str(lumifb)+" fb^{-1}"
        z_ftau.SetTitle("")
	titleText.Draw('sames')
        data_ftau.Draw('sames')
        legendftau.Draw('sames')
        savestr = "ftau.png"
        if doLog == True:
                savestr = "Log"+savestr
        canvasftau.SaveAs(selectiondir+savestr)


make_fakerate(preselisodir, preselnoisodir, preseldir)
make_fakerate(signalisodir, signalnoisodir, signaldir)

