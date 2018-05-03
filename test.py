import bibot



#bibot.plot_country_stats("meta")

#bibot.plot_articles_stats("bibot.log")

#truc = bibot.get_date_from_meta_save("meta/28867810.csv")
#print truc
#bibot.fetch_abstract(28867810)
#bibot.plot_publications_years("meta")



def write_tex_report():
	"""
	IN PROGRESS
	""" 

	## Generate figures
	bibot.plot_publications_years("meta")
	bibot.plot_country_stats("meta")
	bibot.plot_articles_stats("bibot.log")

	## get informations
	keywords_list = []
	total_articles = -1
	selected_articles = 0

	log_file = open("bibot.log", "r")
	cmpt = 0
	for line in log_file:
		line = line.replace("\n", "")
		if(cmpt == 0):
			keywords_list = line.split(";")
		
		elif(cmpt == 1):
			line_in_array = line.split(";")
			total_articles = line_in_array[1]

		elif(line[0] == ">"):

			line_in_array = line.split(";")
			filter_1_info = line_in_array[1].split("=")
			filter_2_info = line_in_array[2].split("=")
			
			if(filter_1_info[1] == "PASSED" and filter_2_info[1] == "PASSED"):
				selected_articles += 1

		cmpt += 1

	log_file.close()

	## create report file
	report_file = open("report.tex", "w")

	## write header
	report_file.write("\\documentclass[a4paper,9pt]{extarticle}\n")
	report_file.write("\\usepackage[utf8]{inputenc}\n")
	report_file.write("\\usepackage[T1]{fontenc}\n")
	report_file.write("\\usepackage{graphicx}\n")
	report_file.write("\\usepackage{xcolor}\n")
	report_file.write("\\usepackage{amsmath,amssymb,textcomp}\n")
	report_file.write("\\everymath{\displaystyle}\n")
	report_file.write("\\usepackage{times}\n")
	report_file.write("\\renewcommand\\familydefault{\sfdefault}\n")
	report_file.write("\\usepackage{tgheros}\n")
	report_file.write("\\usepackage[defaultmono,scale=0.85]{droidmono}\n")
	report_file.write("\\usepackage{multicol}\n")
	report_file.write("\\setlength{\columnseprule}{0pt}\n")
	report_file.write("\\setlength{\columnsep}{20.0pt}\n")
	report_file.write("\\usepackage{geometry}\n")
	report_file.write("\\geometry{\n")
	report_file.write("    a4paper,\n")
	report_file.write("    total={210mm,297mm},\n")
	report_file.write("    left=10mm,right=10mm,top=10mm,bottom=15mm}\n")
	report_file.write("\\linespread{1.3}\n")
	report_file.write("\\makeatletter\n")
	report_file.write("\\renewcommand*{\\maketitle}{%\n")
	report_file.write("\\noindent\n")
	report_file.write("\\begin{minipage}{0.4\\textwidth}\n")
	report_file.write("    \\begin{tikzpicture}\n")
	report_file.write("    \\node[rectangle,rounded corners=6pt,inner sep=10pt,fill=blue!50!black,text width= 0.75\\textwidth] {\\color{white}\\Huge \\@title};\n")
	report_file.write("    \\end{tikzpicture}\n")
	report_file.write("\\end{minipage}\n")
	report_file.write("\\hfill\n")
	report_file.write("\\begin{minipage}{0.55\\textwidth}\n")
	report_file.write("    \\begin{tikzpicture}\n")
	report_file.write("    \\node[rectangle,rounded corners=3pt,inner sep=10pt,draw=blue!50!black,text width= 0.95\\textwidth] {\LARGE \\@author};\n")
	report_file.write("    \\end{tikzpicture}\n")
	report_file.write("\\end{minipage}\n")
	report_file.write("\\bigskip\\bigskip\n")
	report_file.write("}%\n")
	report_file.write("\\makeatother\n")
	report_file.write("\\usepackage[explicit]{titlesec}\n")
	report_file.write("\\newcommand*\\sectionlabel{}\n")
	report_file.write("\\titleformat{\\section}\n")
	report_file.write("{\\gdef\\sectionlabel{}\n")
	report_file.write("\\normalfont\\sffamily\\Large\\bfseries\\scshape}\n")
	report_file.write("{\\gdef\\sectionlabel{\\thesection\\ }}{0pt}\n")
	report_file.write("{\n")
	report_file.write("    \\noindent\n")
	report_file.write("    \\begin{tikzpicture}\n")
	report_file.write("    \\node[rectangle,rounded corners=3pt,inner sep=4pt,fill=blue!50!black,text width= 0.95\\columnwidth] {\\color{white}\\sectionlabel#1};\n")
	report_file.write("    \\end{tikzpicture}\n")		
	report_file.write("}\n")
	report_file.write("\\titlespacing*{\section}{0pt}{15pt}{10pt}\n")	
	report_file.write("\\usepackage{fancyhdr}\n")
	report_file.write("\\makeatletter\n")
	report_file.write("\\pagestyle{fancy}\n")
	report_file.write("\\fancyhead{}\n")
	report_file.write("\\fancyfoot[C]{\\footnotesize \\@date\\ \\ \\@author}\n")	
	report_file.write("\\renewcommand{\\headrulewidth}{0pt}\n")
	report_file.write("\\renewcommand{\\footrulewidth}{0pt}\n")
	report_file.write("\\makeatother\n")
	report_file.write("\\usepackage{tikz}\n")
	report_file.write("\\usetikzlibrary{shapes,arrows}\n")
	report_file.write("\\usepackage{xcolor}\n")
	report_file.write("\\usepackage{chronosys}\n")
	report_file.write("\\usepackage{pgfplots}\n")
	report_file.write("\\usepackage{amsmath,amssymb,textcomp}\n")
	report_file.write("\\everymath{\\displaystyle}\n")

	## write title
	report_file.write("\\title{BIBOT Report}\n")
	report_file.write("\\author{\n")
	report_file.write("\\begin{center}\n")
	report_file.write("   \\vspace{-0.5cm}\n")
	report_file.write("   FIXE DATE\n")
	report_file.write("\\end{center}\n")
	report_file.write("}\n")
	report_file.write("\\date{}\n")

	## write document
	report_file.write("\\begin{document}\n")
	report_file.write("\\maketitle\n")
	report_file.write("\\section{Overview}\n")
	report_file.write("\\subsection{Keywords}\n")

	disposition = "{" + "c" * len(keywords_list) +"}"
	tabular_line = ""
	for keyword in keywords_list:
		tabular_line += keyword + "&"
	tabular_line = tabular_line[:-1]

	report_file.write("\\begin{center}\n")
	report_file.write("\\begin{tabular}"+str(disposition)+"\n")
	report_file.write(str(tabular_line)+"\\\\\n")
	report_file.write("\\end{tabular}\n")
	report_file.write("\\end{center}\n")
	report_file.write("\\subsection{Distributions}\n")
	report_file.write("\\begin{center}\n")
	report_file.write("\\begin{tabular}{cc}\n")
	report_file.write("\\includegraphics[scale=0.5]{images/country_repartition} & \\includegraphics[scale=0.5]{images/years_publications_evolution}\n")
	report_file.write("\\end{tabular}\n")
	report_file.write("\\end{center}\n")
	
	report_file.write("\\section{Filters}\n")
	report_file.write("\\begin{itemize}\n")
	report_file.write("\\item "+str(total_articles)+" articles found\n")
	report_file.write("\\item "+str(selected_articles)+" articles selected\n")
	report_file.write("\\end{itemize}\n")

	report_file.write("\\begin{center}\n")
	report_file.write("\\includegraphics[scale=0.5]{images/plot_stats}\n")
	report_file.write("\\end{center}\n")
	report_file.write("\\end{document}\n")

	## close report file
	report_file.close()


write_tex_report()