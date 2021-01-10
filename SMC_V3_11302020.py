
# This python script was written by a simple biochemist with little formal programming training. Is the code great?
# Nope, but it is indeed functional. The program was designed to accomplish a few quick
# tasks, namely: how many ms2 events occur for a given precursors mass, summing ms2 ions for a given precursor 
# mass, and identifying possible positional isomers along with their respective ratios. As Histone Jim, many  
# of my thought processes started from a TDMS perspective with an analysis of histone proteins in mind. Also, I wrote 
# this code to integrate with software platforms that are freely available  (MASH Explorer, TopPIC, etc.). If
# you'd like to modify this code for your needs, please do so! If you have any questions about it, you can reach
# me at jjp6@stmarys-ca.edu. 

# If using pyinstaller to make windows and mac executable programs, you must downgrade
# matplotlib to version 3.1.0 or else it bugs out.

import numpy as np
import tkinter as tk
import time, sys
import matplotlib
#import matplotlib.pyplot as plt
import matplotlib.backends._tkagg
matplotlib.use("TkAgg")
from matplotlib import style
style.use('ggplot')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from operator import itemgetter, attrgetter
from datetime import datetime
from tkinter import filedialog as fd
from tkinter import messagebox as mb
from tkinter import ttk
from tkinter import *
from collections.abc import Iterable


s_fields = 'Precursor Monoiso. Mass (Da)', 'Start Scan', 'End Scan', 'Mass Tolerance (Da)', 'Histogram Bin Size'
i_fields = 'PTM Isomer Mass (Da)', 'Isomer Mass Tolerance (Da)'
file_types = {'MSALIGN','CSV', 'TDValidator CSV'}
precursor_ions = []
entries = []
i_ents = []
filename = 'temp'

class App(object):

    def __init__(self,master):
        root.geometry("650x700")
        root.wm_title("Search MS1 and Combine MS2 [SMC]")
        self.content = tk.Frame(root, padx=2, pady=2)
        self.content.grid(column=0, row=0, sticky='nsew') 
        self.file_var = StringVar(root)
        self.file_var.set('Select Type/File') 
        self.popupMenu = OptionMenu(self.content, self.file_var, *file_types)
        Label(self.content, text="Select data type and open file:").grid(row = 0, column = 0)
        self.popupMenu.grid(row = 0, column =1)
        self.file_var.trace('w',self.makeform)

        self.test_isomers = IntVar()
        self.isomer_btn = tk.Checkbutton(self.content, text="Isomer Search (Optional)", variable=self.test_isomers)

        self.b1 = tk.Button(self.content, text='Run',  fg="green", command=(lambda : self.fetch()))
        self.b2 = tk.Button(self.content, text='Quit', command=root.destroy)
        self.b2.grid(row=6,column=0)
        self.e=StringVar()
        self.e.set("                                                 ")
        self.loading = tk.Label(self.content, textvariable=self.e)
        self.loading.grid(row=8, columnspan=4, sticky="ew")
        
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        self.content.columnconfigure(0, weight=2)
        self.content.columnconfigure(1, weight=2)  

       
    def clear(self): 
        #self.file_var = StringVar(root)
        #self.file_var.set('Select Type/File') 
        #self.popupMenu = OptionMenu(self.content, self.file_var, *file_types)
        Label(self.content, text="Select data type and open file:").grid(row = 0, column = 0)
        self.popupMenu.grid(row = 0, column =1)
        #self.file_var.trace('w',self.makeform)

        self.test_isomers = IntVar()
        self.isomer_btn = tk.Checkbutton(self.content, text="Isomer Search (Optional)", variable=self.test_isomers)

        #self.b1 = tk.Button(self.content, text='Run',  fg="green", command=(lambda : self.fetch()))
        #self.b2 = tk.Button(self.content, text='Quit', command=root.destroy)
        self.b2.grid(row=6,column=0)
        self.e=StringVar()
        self.e.set("                                                 ")
        self.loading = tk.Label(self.content, textvariable=self.e)
        self.loading.grid(row=8, columnspan=4, sticky="ew")
        
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        self.content.columnconfigure(0, weight=2)
        self.content.columnconfigure(1, weight=2)
        
    
    def update_progress(self,progress):
        barLength = 15 # Modify this to change the length of the progress bar
        status = ""
        if isinstance(progress, int):
            progress = float(progress)
        if not isinstance(progress, float):
            progress = 0
            status = "error: progress var must be float\r\n"
        if progress < 0:
            progress = 0
            status = "Halt...\r\n"
        if progress >= 1:
            progress = 1
            status = "Done...\r\n"
        block = int(round(barLength*progress))
        progresstext = "\rLoading (%): [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), int(progress*100), status)
        self.e.set(progresstext)
        #sys.stdout.write(progresstext)
    
    def display_output(self,filenames):
        outputframe = tk.Tk()
        outputframe.geometry("800x500")
        outputframe.wm_title("SMC Results")
        j=0

        for i in filenames:
            with open(i, "r") as f:
                data = f.read()
                data = data.replace("{}", "")
                l = tk.Label(outputframe,
                             text="Data from file name:\n"+i,
                             font = 'Times 12',
                             fg = 'green',
                             bd = 1,
                             relief = 'solid')
                t = tk.Text(outputframe,
                                 font = 'Times 10')
                t.insert(tk.END, data)
                t.grid(column=j,row=1,padx=3,sticky='news')
                l.grid(column=j,row=0)
                scrollb = ttk.Scrollbar(outputframe, command=t.yview)
                scrollb.grid(row=1, column=j+1, sticky='news')
                t['yscrollcommand'] = scrollb.set
                outputframe.columnconfigure(j, weight=1)
                outputframe.columnconfigure(j+1,minsize=50, weight=1)
                outputframe.rowconfigure(j+1, weight=1)
            j+=2


    #remove all characters from line function
    def rmv_chars(self,string):
        getVals = list([val for val in string 
            if (val.isnumeric() or val==".")]) 

        return "".join(getVals)


    def combine_ms2(self,array, mass, scanfirst, scanlast, tolerance):
        sent_array = np.array(array, dtype=object)
        new_array = []
        temp = []
        temp = np.array(temp,dtype=object)
        sorted_array = []
        sorted_array= np.array(sorted_array, dtype=object)

        #Find and store the number of activation events
        for i in range(len(array)):
            temp = np.append(temp, array[i,3])
        activation_types = np.unique(temp)

        if (len(sent_array)) == 0:
            return mb.showerror("No precursor ions within",tolerance,"Da of",mass,"to combine.")

        #filter based off search criteria
        for i in range(len(sent_array)):
            if (sent_array[i][6] < (mass+tolerance)) and (sent_array[i][6] > (mass-tolerance)) and (sent_array[i][1] >= scanfirst) and (sent_array[i][1] <= scanlast):
                new_array.extend((sent_array[i][3],sent_array[i][6],sent_array[i][7],sent_array[i][1]))
        new_array = np.reshape(new_array,(int(len(new_array)/4),4))
        #new_array.dtype = object

        #sort by activation
        #final sort array is organized as [activation, precursor mass, [ms2 info],scan number]
        for i in activation_types:
            for j in range(len(new_array)):
                if i == new_array[j,0]:
                    sorted_array = np.append(sorted_array,new_array[j])
        sorted_array = np.reshape(sorted_array,(int(len(sorted_array)/4),4))

        return sorted_array


    #This function takes an activation type-sorted array of [Activation type, Prec Mass, [MS2 ion info],Scan number] 
    #and generates an output file per activation type. This file includes the searched mass, scans that were summed,
    #and the mass of each fragment ion from MS2 data from the searched, parameterized precursor mass.
    def output_file(self,sorted_array, mass):
        ms2_ions = []
        temp_ions = []
        date = datetime.today().strftime('%Y%m%d_%H%M%S')
        filenames = []
        a = []

        if (len(sorted_array)) == 0:
            mb.showerror("Warning", "Based on your search criteria, no results were generated.")
        else:
            activation = str(sorted_array[0][0])
            while len(sorted_array) != 0:
                ms2_ions.append(activation)
                filename = str(int(mass))+"_"+(activation).replace('ACTIVATION=','')+"_"+date+".txt"
                with open(filename,"a") as out_file:
                    filenames.append(filename)    
                    out_file.write(activation+"\nSearched Mass="+str(mass)+"\nScans summed: ")
                    for i in sorted_array:
                        if str(i[0]) == activation:
                            out_file.write(str(i[3])+", ")
                    out_file.write(" \nPrecursor ions selected:")
                    for i in sorted_array:
                        if str(i[0]) == activation:
                            out_file.write(str(i[1])+", ")
                    out_file.write(" \n")
                    for i in sorted_array:
                        if str(sorted_array[0][0]) != activation:
                            activation = str(sorted_array[0][0])
                            break
                        if i[0] == activation:
                            if isinstance(i[2],Iterable):
                                for j in i[2]:
                                    out_file.write("\n"+str(j[0]))
                                    temp_ions.append(j)
                            a.append(temp_ions)
                            sorted_array = np.delete(sorted_array,0,0)
                        temp_ions = []
                    ms2_ions.append(a)
                    a = []

            #Do isomer search. ms2_ions has the following structure [((activation),([([a,b,c],)],),] where
            # a,b,c is the monoisotopic mass, intensity and charge of an ion within a single scan. 
            # ms2_ions[1] is therefore a collection of all returned scans for a searched precursor mass from an
            # activation type (in ms2_ions[0]). For example, ms2_ions[1][0] would
            # have all the the ions for a single scan and ms2_ions[1][0][0] would have information on just
            # one ion from one scan from one activation type.
            ms2_ions = np.array(ms2_ions, dtype=object)
            if self.test_isomers.get() == 1:
                self.isomer_search(ms2_ions,mass)        

            self.display_output(filenames)       
            mb.showinfo("Run Completed", "The run has completed and generated an output file(s).")

            
    def makegraph(self,precursors, hist_bin_size):
        font = {'family' : 'DejaVu Sans','weight' : 'bold','size'   : 12}

        bin_size = int(hist_bin_size.get())

        matplotlib.rc('font', **font)
        f = Figure(figsize=(1,4), dpi=100, tight_layout=True)
        a = f.add_subplot()
        a.hist(precursors, bin_size, density=False, facecolor='g', alpha=0.75)

        hist_bin_size.bind('<Return>', (lambda event: self.makegraph(precursors, hist_bin_size)))

        a.set_xlabel('Precursor Mass (Da)')
        a.set_ylabel('Number of MS2 Events')
        a.set_title('Histogram of Precursor Masses')

        canvas = FigureCanvasTkAgg(f, self.content)
        canvas.draw()
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.grid(row=8, columnspan=4, sticky='news')

        toolbar_frame = tk.Frame(self.content) 
        toolbar_frame.grid(row=7, columnspan=4) 
        toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
            
        return


    def fetch(self):
        data = []
        data_isomers = []
        filetype = self.file_var.get()
        error = False

        if filetype == 'MSALIGN': 
            for entry in entries: 
                if entry[1].get() == "":
                    error = self.callback("fail", entry[0])
                elif float(entry[1].get()) > 100 and entry[0] == "Mass Tolerance (Da)":
                    mb.showerror("Warning", "Mass Tolerance has to be less than 100 Da.")
                    error = True

        if filetype == 'CSV':
            for entry in i_ents: 
                if entry[1].get() == "":
                    error = self.callback("fail", entry[0])
            if error:
                return
            else:
                self.process_CSV()
        
        if error:
            return

        else:         
            for entry in entries:
                data.append(float(entry[1].get()))

            if filetype == "MSALIGN":
                self.callback("data", data)


    def callback(self,param, data):
        self.clear()
        if param == "MSALIGN":
            self.b1.grid(row=6,column=1)
            scan = []
            entries[3][1].insert(0,1)
            entries[4][1].insert(0,1000)

            self.makegraph(self.process(data), entries[4][1])

            #Autopopulate scan range (based on loaded MSALIGN file) and mass tolerance    
            for scans in precursor_ions:
                scan.append(scans[1])
            entries[1][1].insert(0,np.amin(scan))
            entries[2][1].insert(0,np.amax(scan))

        if param == "CSV":
            self.b1.grid(row=6,column=1)
            self.test_isomers.set(1)

        if param == "fail":
            mb.showerror("Warning", "Error: the "+str(data)+" field is empty.")
            return True

        if param == "data":
            sorted_array = self.combine_ms2(precursor_ions, data[0], data[1], data[2], data[3])
            self.output_file(sorted_array, data[0])


    def makeform(self,*args):
        global entries
        global i_ents
        global filename
        global precursor_ions
        
        precursor_ions = []
        entries = []
        i_ents = []
        i=1
        
        filetype = self.file_var.get()
        
        if filetype == 'MSALIGN':
            extension = '.msalign'
        else:
            extension = '.csv'
            
        self.isomer_btn.grid(row=0,column=2, columnspan=2, sticky='news')
        
        filename = fd.askopenfilename(filetypes=[(filetype,extension),('All files','*.*')],
                                      title='Please select the '+filetype+' MS data file')


        if filetype == 'MSALIGN':
            for field in s_fields:
                lab = tk.Label(self.content, width=20, text=field, anchor='w')
                ent = tk.Entry(self.content, width=8)
                lab.grid(row=i, column=0)
                ent.grid(row=i, column=1)
                entries.append((field, ent))
                i+=1
            i=1

        for field in i_fields:
            lab = tk.Label(self.content, width=20, text=field, anchor='w')
            ent = tk.Entry(self.content, width=8)
            lab.grid(row=i, column=2)
            ent.grid(row=i, column=3)
            i_ents.append((field, ent))
            i+=1          
        
        if filetype == 'TDValidator CSV':
            search_params = fd.askopenfilename(filetypes=[('Search parameter','.csv')],title='Please select the isomer criteria file')
            if search_params == "":
                self.makeform
            else:
                self.tdvalidator_search(search_params)
        
        if filename == "":
            print('no file selected')
            self.makeform
        else:
            self.callback(filetype, filename)


    def process(self,name):
        global precursor_ions
        precursor_ions = []
        msfile = name
        temp_array = [0,0,0,0,0,0,0,0]
        precursors = []
        ms2_ions = []
        ms2_ions = np.array(ms2_ions)
        ms2events = 0
        counter = 0
        convert = []
        lines = []
        
        with open(msfile) as fp:
            data = fp.read()

        for i in data:
            if i == ',' or i == '\n':
                link = "".join(convert)
                link.strip()
                if len(link) > 1:
                    lines.append(link)
                    convert = []
            else:
                convert.append(i)
        progress_bar_lines = len(lines)
        
        
        while (len(lines)-1)>0:
            text = str(lines[0])
            if text.startswith('ID='):
                temp_array[0] = int(self.rmv_chars(text))
                counter+=1
                self.update_progress((1-(len(lines)/progress_bar_lines))) #do progress update here b/c of time
                self.loading.configure(textvariable=self.e)
                root.update()
            elif text.startswith('SCANS='):
                temp_array[1] = int(self.rmv_chars(text))
                counter+=1
            elif text.startswith('RETENTION_TIME='):
                temp_array[2] = float(self.rmv_chars(text))
                counter+=1
            elif text.startswith('ACTIVATION='):
                temp_array[3] = str(text)
                counter+=1
            elif text.startswith('PRECURSOR_MZ='):
                temp_array[4] = float(self.rmv_chars(text))
                counter+=1
            elif text.startswith('PRECURSOR_CHARGE='):
                temp_array[5] = int(self.rmv_chars(text))
                counter+=1
            elif text.startswith('PRECURSOR_MASS='):
                temp_array[6] = float(self.rmv_chars(text))
                counter+=1  
            elif text[0].isdigit():
                counter+=1
                while lines[0]!='END IONS':
                    ms2_ions = np.append(ms2_ions,[float(s) for s in lines[0].split("\t")])
                    del lines[0]
                ms2_ions = np.reshape(ms2_ions,(int(len(ms2_ions)/3),3))
                temp_array[7] = ms2_ions
                precursor_ions.append(temp_array)
                temp_array = [0,0,0,0,0,0,0,0]
                ms2_ions = []
                ms2events+=1
            del lines[0]

        self.loading.grid_forget()
        
        #This is very important: Here is the final format of the np array:
        #For each precursor ion, the format is [ID,SCANS,RETENTION_TIME,ACTIVATION,PRECURSOR_MZ,PRECURSOR_CHARGE,PRECURSOR_MASS,[ARRAY OF MS2 IONS]]
        #The array of MS2 ions is [Mass,Intensity,Charge]
        precursor_ions = np.array(precursor_ions, dtype=object)

        #Autopopulate Precursor Mass if the MSALIGN file has only a single precursor ion (i.e., Thrash deconv from Mash Explorer)
        if ms2events == 1:
            entries[0][1].insert(0,precursor_ions[0][6])

        for i in precursor_ions:
            precursors.append(i[6])

        precursors = np.array(precursors)
        root.wm_title(name)
        return(precursors)

    def process_CSV(self):
        mass = 0
        activation = 'ACTIVATION=CSV'
        msfile = filename
        convert = []
        pool = []
        scan_ions = []
        organized_CSV = []
        j=0

        #input eThrash or csv file and convert it to string list
        with open(msfile) as fp:
            data = fp.read()

        for i in data:
            if i == ',' or i == '\n':
                link = "".join(convert)
                pool.append(link)
                convert = []
            else:
                convert.append(i)
        pool = np.array(pool, dtype=object)
        pool = np.delete(pool,0,0)

        while len(pool) > 0:
            if pool[j] == "":
                pool = np.delete(pool,0,0)
            elif pool[j][0].isnumeric():
                a = [float(pool[j]),float(pool[j+1]),int(pool[j+2])]
                scan_ions.append(a)
                pool = np.delete(pool,[0,1,2],None)
            else:
                pool = np.delete(pool,0,0)
                organized_CSV.append(scan_ions)
                scan_ions = []
            a=[]

        organized_CSV = np.array((activation,organized_CSV), dtype=object)
        self.isomer_search(organized_CSV,mass)

    def tdvalidator_search(self,search_params):
        pool = []
        convert = []
        matched_fragments = []
        date = datetime.today().strftime('%Y%m%d_%H%M%S')
        chargestates = [1,2,3]
        mz_search_list = []
        match_list = []
        proton = 1.007825
        tolerance_ppm = 20
        
        
        with open(search_params) as pfp:
            data_params = pfp.read()
            
        with open(filename) as dfp:
            data_tdvalidator = np.loadtxt(dfp, dtype=float, delimiter=",", unpack=False)

        for idx, val in enumerate(data_params):
            if val == ',' or val == '\n' or idx == (len(data_params)-1):
                link = "".join(convert)
                pool.append(link)
                convert = []
            else:
                convert.append(val)
        pool = np.array(pool, dtype=object)
        
        #reshaped parameters is now [(amino acid, position, acetylstate, mass)]
        parameters = np.reshape(pool,(int(len(pool)/4),4))
        
        #convert fragment monoisotopic masses to m/z values based on chosen charge states
        for i in parameters:
            temp = []
            temp.extend((i[0],i[1],i[2]))
            for j in chargestates:
                temp.append((float(i[3])+j*proton)/(j))
            mz_search_list.append(temp)

        for j in range(len(chargestates)):
            for i in data_tdvalidator:
                k=0
                observed = i[0]
                while k < len(mz_search_list):
                    temp_list = []
                    theoretical = mz_search_list[k][j+3]
                    delta = abs(theoretical - observed)
                    ppm_error = delta/theoretical * (10 ** 6)
                    if ppm_error <= tolerance_ppm:
                        temp_list.extend((mz_search_list[k][0],int(mz_search_list[k][1]),mz_search_list[k][2],
                                         theoretical, observed,'z= '+str(chargestates[j]), i[1]))
                        match_list.append(temp_list)
                    k+=1
        match_list = sorted(match_list, key=itemgetter(1,2))
        
        tdval_filename = 'TDvalidator_out_'+date+"_search.csv"
        with open(tdval_filename,"a") as out_file:
            out_file.write("Residue, Position, PTM, Theoretical Mass, Observed Mass, Charge state, Intensity")
            for i in match_list:
                line = str(i)[1:-1]
                out_file.write("\n"+line)
        
    def isomer_search(self,ms2_ions,mass):
        date = datetime.today().strftime('%Y%m%d_%H%M%S')
        matched_fragments = []
        massdiff, tolerance = float(i_ents[0][1].get()),float(i_ents[1][1].get())

        for activation_type in range(0,len(ms2_ions),2):
            filename = str(int(mass))+"_"+str(ms2_ions[activation_type]).replace('ACTIVATION=','')+"_"+date+"_isomers.txt"
            with open(filename,"a") as out_file:
                for scan in range(len(ms2_ions[activation_type+1])):
                    while len(ms2_ions[activation_type+1][scan]) > 0:
                        for i in range(len(ms2_ions[activation_type+1][scan])-1):
                            fmass = float(ms2_ions[activation_type+1][scan][0][0])
                            nmass = float(ms2_ions[activation_type+1][scan][i][0])
                            fint = float(ms2_ions[activation_type+1][scan][0][1])
                            nint = float(ms2_ions[activation_type+1][scan][i][1])
                            fcstate = int(ms2_ions[activation_type+1][scan][0][2]) 
                            ncstate = int(ms2_ions[activation_type+1][scan][i][2])

                            if fmass >= nmass:
                                test = fmass - nmass
                            else:
                                test = nmass - fmass
                            if (massdiff - tolerance) <= test <= (massdiff + tolerance) and fcstate == ncstate:
                                if fmass >= nmass:
                                    ratio = (fint/(fint+nint)) * 100
                                    a = [fmass,nmass,fcstate,test,ratio]
                                else:
                                    ratio = (nint/(fint+nint)) * 100
                                    a = [nmass,fmass,fcstate,test,ratio]
                                matched_fragments.append(a)
                        ms2_ions[activation_type+1][scan] = np.delete(ms2_ions[activation_type+1][scan],0,0)

                matched_fragments = sorted(matched_fragments, key=lambda x: x[0])  
                for i in matched_fragments:
                    out_file.write("\n"+str(i))
                matched_fragments = []

        mb.showinfo("Run Completed", "The run has completed and generated an output file(s).")


root = tk.Tk()
app = App(root)        
root.mainloop()
