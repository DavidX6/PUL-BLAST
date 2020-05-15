#! /usr/bin/env python
#  -*- coding: utf-8 -*-
#
# GUI module generated by PAGE version 5.0.3
#  in conjunction with Tcl version 8.6
#    Apr 04, 2020 05:55:16 PM CEST  platform: Windows NT

from tkinter import scrolledtext
from tkinter import filedialog
import webbrowser
import os

import tkinter as tk
import tkinter.ttk as ttk

py3 = True

import gui_support
from main import *


def vp_start_gui():
    '''Starting point when module is the main routine.'''
    global val, w, root
    root = tk.Tk()
    gui_support.set_Tk_var()
    top = Toplevel1(root)
    gui_support.init(root, top)
    root.mainloop()


w = None


def create_Toplevel1(rt, *args, **kwargs):
    '''Starting point when module is imported by another module.
       Correct form of call: 'create_Toplevel1(root, *args, **kwargs)' .'''
    global w, w_win, root
    # rt = root
    root = rt
    w = tk.Toplevel(root)
    gui_support.set_Tk_var()
    top = Toplevel1(w)
    gui_support.init(w, top, *args, **kwargs)
    return (w, top)


def destroy_Toplevel1():
    global w
    w.destroy()
    w = None


class Toplevel1:
    global genomeResults
    genomeResults = {}
    global modularityDict
    modularityDict = {}

    def fillModularity(self):
        dct = {}
        f = open("modularityDescription.txt", "r", encoding="utf-8")
        currentKey = ""
        for line in f:
            if line[0].isdigit(): currentKey = line.replace("\n", "")
            else: dct[currentKey] = line.replace("\n", "").replace(u'\xa0', ' ')
        return dct

    def validUserEval(self, userEval):
        try: return float(userEval)
        except: return 0.001

    def validUserOut(self, userOut):
        try:
            a = int(userOut)
            if a >= 0 and a <= 11:
                return a
            else: return 5
        except: return 5

    def validUserDist(self, userDist):
        try:
            a = int(userDist)
            if a >= 1:
                return a
            else: return 10
        except: return 10

    def browseButtonSubmit(self, type):
        userSequence = self.Entry1.get()
        self.Entry1.delete(0, tk.END)
        if len(userSequence) < 1 or userSequence == "processing" or userSequence == "No sequence":
            self.Entry1.insert(0, "No sequence")
        else:
            self.Entry1.insert(0, "processing")
            userOut = self.validUserOut(self.EntryOutfmt.get())
            userEval = self.validUserEval(self.EntryEval.get())
            userDist = self.validUserDist(self.EntryDist.get())
            print(userOut, userEval, userOut)
            self.ComboBox["values"] = ["None"]
            self.ComboBox.current(0)
            self.Text1.delete(1.0, tk.END)
            self.Text2.delete(1.0, tk.END)
            self.modularityDict = self.fillModularity()
            if type not in ["gbk", "fasta"]:
                self.Entry1.delete(0, tk.END)
                self.Entry1.insert(0, "no type selected")
            else:
                if type == "fasta": genomeSearch(userSequence)
                else: gbkGenomeSearch(userSequence)
                proteinBLAST("queryTempSequence.txt", eval=userEval)
                resList = searchPULs(maxDist=userDist)
                self.genomeResults = resList
                self.Entry1.delete(0, tk.END)
                self.Entry1.insert(0, "processed")
                self.fillOptions()
                if userOut != 5: proteinBLAST("queryTempSequence.txt", eval=userEval, format=userOut,
                                              outfile="UserRequestedBLASTres.txt")

    def writePULs(self, substrate):
        self.Text2.delete(1.0, tk.END)
        for key in self.modularityDict.keys():
            if substrate == key[key.index(" ")+1 : ]:
                self.Text2.insert("insert", key + " " + self.modularityDict[key] + "\n")

        self.Text1.delete(1.0, tk.END)
        self.Text1.insert("insert", str(len(self.genomeResults)) + " candidate PULs found" + "\n")
        for pul in self.genomeResults:
            self.Text1.insert("insert", "-PUL match details-" + "\n")
            for record in pul:
                cnt = 0
                for alignment in record.alignments:
                    currentSubstrate = alignment.hit_def.split("|")[3].strip()
                    if currentSubstrate != substrate and substrate != "Everything": break
                    if cnt == 0:
                        self.Text1.insert("insert", record.query + "\n")
                        cnt += 1
                    self.Text1.insert("insert", "\t" + alignment.hit_def + "\n")
                    for hsp in alignment.hsps:
                        self.Text1.insert("insert",
                                          "\t" + "Bit score: " + str(hsp.bits) + ", evalue: " + str(hsp.expect) + "\n")
                        self.Text1.insert("insert", "\t" + "Identities: " +
                                          str("{:.2f}".format(hsp.identities*100/hsp.align_length)) + "% (" +str(hsp.identities) + "), gaps: " +
                                          str(hsp.gaps) +"\n")
                        self.Text1.insert("insert", "\t" + "Query range: " + str(hsp.query_start) + "-" + str(
                            hsp.query_end) + "\n")
                        self.Text1.insert("insert", "\t" + "Match range: " + str(hsp.sbjct_start) + "-" + str(
                            hsp.sbjct_end) + "\n")
                        self.Text1.insert("insert", "\t" + hsp.query + "\n")
                        self.Text1.insert("insert", "\t" + hsp.match + "\n")
                        self.Text1.insert("insert", "\t" + hsp.sbjct + "\n")
                        self.Text1.insert("insert", "\n")
            self.Text1.insert("insert", "\n")

    def fillOptions(self):
        substrates = set()
        for pul in self.genomeResults:
            for record in pul:
                for alignment in record.alignments: substrates.add(alignment.hit_def.split("|")[3].strip())
        temp = []
        for string in substrates:
            if string.strip() == "General": continue
            temp.append(string)
        temp.append("Everything")
        self.ComboBox["values"] = temp

    def showSubstrate(self, *args):
        val = self.ComboBox["values"][self.ComboBox.current()]
        if val.lower() != "none": self.writePULs(val)

    def askopenfile(self):
        a = tk.filedialog.askopenfilename(title="Select file")
        self.Entry1.delete(0, tk.END)
        self.Entry1.insert(0, a)
        gui_support.selectedButton.set(0)
        gui_support.selectedGenus.set(0)

    def extraFeatures(self, command):
        if command == "export":
            f = tk.filedialog.asksaveasfile(mode='w', defaultextension=".txt")
            if f is not None:
                results = self.Text1.get("1.0", tk.END)
                f.write(results)
                f.close()
        elif command == "open": webbrowser.open('file://' + os.path.realpath("javascript/index.html"), new=1)

    def __init__(self, top=None):
        '''This class configures and populates the toplevel window.
           top is the toplevel containing window.'''

        top.geometry("893x659")
        top.minsize(148, 1)
        top.maxsize(1924, 1055)
        top.resizable(1, 1)
        top.title("PUL BLAST")
        # file selection
        self.Label1 = ttk.Label(top)
        self.Label1.place(relx=0.032, rely=0.029, height=37, width=220)
        self.Label1.configure(text="Select sequence file or enter path")

        self.Button1 = ttk.Button(top, command=self.askopenfile)
        self.Button1.place(relx=0.82, rely=0.09, height=33, width=80)
        self.Button1.configure(text="Browse")

        self.Entry1 = ttk.Entry(top)
        self.Entry1.place(relx=0.034, rely=0.09, height=34, relwidth=0.761)

        # BLAST parameters
        self.EntryOutfmt = ttk.Entry(top)
        self.EntryOutfmt.place(relx=0.12, rely=0.155, height=20, relwidth=0.1)
        self.Label2 = ttk.Label(top)
        self.Label2.place(relx=0.034, rely=0.15, height=30, width=50)
        self.Label2.configure(text="Outfmt*:")

        self.EntryEval = ttk.Entry(top)
        self.EntryEval.place(relx=0.32, rely=0.155, height=20, relwidth=0.1)
        self.Label3 = ttk.Label(top)
        self.Label3.place(relx=0.234, rely=0.15, height=30, width=50)
        self.Label3.configure(text="Evalue*:")

        self.EntryDist = ttk.Entry(top)
        self.EntryDist.place(relx=0.57, rely=0.155, height=20, relwidth=0.1)
        self.Label5 = ttk.Label(top)
        self.Label5.place(relx=0.434, rely=0.15, height=30, width=80)
        self.Label5.configure(text="Max distance*:")

        self.Label4 = ttk.Label(top)
        self.Label4.place(relx=0.811, rely=0.15, height=30, width=60)
        self.Label4.configure(text="(*optional)")

        # input type button
        self.ButtonF = ttk.Button(top, command=lambda: self.browseButtonSubmit("fasta"))
        self.ButtonF.place(relx=0.25, rely=0.250, height=40, width = 150)
        self.ButtonF.configure(text="FASTA genome")

        self.ButtonG = ttk.Button(top, command=lambda: self.browseButtonSubmit("gbk"))
        self.ButtonG.place(relx=0.5, rely=0.250, height=40, width = 150)
        self.ButtonG.configure(text="GenBank genome")

        # results display
        self.substratesList = ["None"]
        self.Label4 = ttk.Label(top)
        self.Label4.place(relx=0.29, rely=0.33, height=40, width = 150)
        self.Label4.configure(text="Select substrate:")

        self.ComboBox = ttk.Combobox(top, state="readonly")
        self.ComboBox["values"] = self.substratesList
        self.ComboBox.current(0)
        self.ComboBox.bind("<<ComboboxSelected>>", self.showSubstrate)
        self.ComboBox.place(relx=0.5, rely=0.330, height=40, width=150)

        self.Text1 = scrolledtext.ScrolledText(top, wrap="none")
        self.Text1.place(relx=0.034, rely=0.425, relheight=0.45, relwidth=0.750)
        self.Text1.configure(background="white")
        self.Text1.configure(font="Consolas 10")
        self.textHsb = ttk.Scrollbar(self.Text1, orient="horizontal", command=self.Text1.xview)
        self.textHsb.pack(side="bottom", fill="x")
        self.Text1.configure(xscrollcommand=self.textHsb.set)

        self.Text2 = tk.scrolledtext.ScrolledText(top, wrap="none")
        self.Text2.place(relx=0.034, rely=0.885, relheight=0.08, relwidth=0.750)
        self.Text2.configure(background="white")
        self.Text2.configure(font="Consolas 10")
        self.textHsb1 = ttk.Scrollbar(self.Text2, orient="horizontal", command=self.Text2.xview)
        self.textHsb1.pack(side="bottom", fill="x")
        self.Text2.configure(xscrollcommand=self.textHsb1.set)

        # extra features
        self.ButtonExport = ttk.Button(top, command=lambda: self.extraFeatures("export"))
        self.ButtonExport.place(relx=0.82, rely=0.425, height=40, width=100)
        self.ButtonExport.configure(text="Export")
        self.ButtonOpen = ttk.Button(top, command=lambda: self.extraFeatures("open"))
        self.ButtonOpen.place(relx=0.82, rely=0.50, height=40, width=100)
        self.ButtonOpen.configure(text="Browse BLAST \n       results")


if __name__ == '__main__':
    vp_start_gui()
