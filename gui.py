#coding: utf8 


from kivy.app import App
from kivy.uix.pagelayout import PageLayout

from kivy.uix.filechooser import FileChooser
from kivy.uix.filechooser import FileChooserIconLayout
from kivy.uix.filechooser import FileChooserListLayout
from kivy.uix.filechooser import FileChooserIconView

from kivy.properties import NumericProperty, ObjectProperty
from kivy.uix.screenmanager import ScreenManager, Screen
from os.path import sep, expanduser, isdir, dirname

from kivy.uix.boxlayout import BoxLayout
from kivy.lang import Builder

import os
import unidecode
import glob
import time
from time import gmtime, strftime
import datetime
import shutil
import threading
import random
import bibot



Builder.load_string("""
<Screen1App>:
    myslider: slider
    user_request: my_user_input
    PageLayout:
        
        BoxLayout:
            canvas:
                Color:
                    rgb: 1, 1, 1, 1
                Rectangle:
                    pos: self.pos
                    size: self.size

            orientation: "vertical"
            Image:
                source: 'images/image.jpg'
                keep_ratio: True
                allow_stretch: True
                opacity: 1
                pos_hint: {'center_x': 0.5, 'center_y': 1}

            TextInput:
                id: my_user_input
                text: "machine learning;SjS;big data"
                size_hint_y: .2
                size_hint_x: .4
                pos_hint: {'center_x': 0.5, 'center_y': 0.5}

            BoxLayout:
                spacing: 150
                padding: 40
                orientation: "horizontal"
                Button:
                    text: "Confirm Request"
                    size_hint_y: .4
                    on_press: root.manager.current = 'Screen2'
                    on_press: root.test_user_request()
                    on_press: root.run_alpha()


        BoxLayout:
            orientation: "vertical"
            canvas:
                Color:
                    rgba: .1, .1, .1, 1
                Rectangle:
                    pos: self.pos
                    size: self.size
            Label:
                markup: True
                text: "Option Frame"
                color: 0, 0, 0, 1
                outline_color: 0, 0.5, 0.5, 1
                font_size: 30
            GridLayout:
                spacing: 50
                padding: 10
                cols: 2
                rows: 2
                Button:
                    text: "bouton 1"
                    on_press: root.manager.current = 'Screen2'
                FileChooser:
                    id: fc
                    FileChooserIconLayout
                    FileChooserListLayout
                
                Slider:
                    id: slider
                    min: 0
                    max: 100
                    value: 25
                
                GridLayout:
                    spacing: 50
                    padding: 10
                    cols: 2
                    rows: 1
                
                    CheckBox:
                        active: True
                    Switch:
                        active: True

            AnchorLayout:
                anchor_x: 'center'
                anchor_y: 'bottom'
                padding: 50
                Button:
                    size_hint: 0.3, 0.4
                    text: "Valider"


        BoxLayout:
            spacing: 25
            padding: 20
            orientation: "vertical"
            canvas:
                Color:
                    rgba: .2, .2, .2, 1
                Rectangle:
                    pos: self.pos
                    size: self.size
            FileChooserListView:
                id: filechooser
                on_selection: root.selected(filechooser.selection)
            Button
                size_hint: 0.3, 0.4
                text: "Load configuration file"
                on_release: root.open(filechooser.path, filechooser.selection)

<Screen2App>:
	progress_bar: pb
	result_button: go_button
    run_button: run_button
    status_message: status
    joke_message: joke
    BoxLayout:
        canvas:
            Color:
                rgb: 1, 1, 1, 1
            Rectangle:
                pos: self.pos
                size: self.size

        orientation: "vertical"
        
        
        Label:
            id: status
            markup: True
            text: "Status"
            color: 0, 0, 0, 1
            outline_color: 0, 0.5, 0.5, 1
            font_size: 20

        AnchorLayout:
            anchor_x: 'center'
            anchor_y: 'bottom'
            padding: 30
            ProgressBar:
                id: pb
                size_hint_x: .5
                size_hint_y: None
                height: '48dp'
                value: 0
        Label:
            id: joke
            markup: True
            color: 0, 0, 0, 1
            outline_color: 0, 0.5, 0.5, 1
            font_size: 15

        AnchorLayout:
            anchor_x: 'center'
            anchor_y: 'bottom'
            padding: 50
            Button:
                id: run_button
                size_hint: 0.5, 2.5
                text: "Run"
                on_release: root.run_articles_selection_wrapper()
    		
            Button:
    			id: go_button
    			size_hint_x: 0
    			size_hint_y: 0
    			opacity:0
    			disabled: True
    			on_press: root.manager.current = 'Screen3'

<Screen3App>:

	PageLayout:

	    BoxLayout:
	        canvas:
	            Color:
	                rgb: 1, 1, 1, 1
	            Rectangle:
	                pos: self.pos
	                size: self.size

	        orientation: "vertical"
	        Label:
	            markup: True
	            text: "Screen 3 - Results summary"
	            color: 0, 0, 0, 1
	            outline_color: 0, 0.5, 0.5, 1
	            font_size: 30
            Image:
                source: 'images/years_publications_evolution.png'
                keep_ratio: True
                allow_stretch: True
                opacity: 1
                pos_hint: {'center_x': 0.5, 'center_y': 1}
            BoxLayout:
                canvas:
                    Color:
                        rgb: 1, 1, 1, 1
                    Rectangle:
                        pos: self.pos
                        size: self.size

                orientation: "horizontal"
                Image:
                    source: 'images/years_publications_evolution.png'
                    keep_ratio: True
                    allow_stretch: True
                    opacity: 1
                    pos_hint: {'center_x': 0.5, 'center_y': 1}
                Image:
                    source: 'images/country_repartition.png'
                    keep_ratio: True
                    allow_stretch: True
                    opacity: 1
                    pos_hint: {'center_x': 0.5, 'center_y': 1}

	    BoxLayout:
	        canvas:
	            Color:
	                rgb: 1, 1, 1, 1
	            Rectangle:
	                pos: self.pos
	                size: self.size

	        orientation: "vertical"
	        Label:
	            markup: True
	            text: "Country Repartition"
	            color: 0, 0, 0, 1
	            outline_color: 0, 0.5, 0.5, 1
	            font_size: 30
            Image:
                source: 'images/country_repartition.png'
                keep_ratio: True
                allow_stretch: True
                opacity: 1
                pos_hint: {'center_x': 0.5, 'center_y': 1}

	    BoxLayout:
	        canvas:
	            Color:
	                rgb: 1, 1, 1, 1
	            Rectangle:
	                pos: self.pos
	                size: self.size

	        orientation: "vertical"
	        Label:
	            markup: True
	            text: "Screen 3 - Results 2"
	            color: 0, 0, 0, 1
	            outline_color: 0, 0.5, 0.5, 1
	            font_size: 30
        


<MyWidget>:
    id: my_widget
    Button
        text: "open"
        on_release: my_widget.open(filechooser.path, filechooser.selection)
    FileChooserListView:
        id: filechooser
        on_selection: my_widget.selected(filechooser.selection)
""")

class MyWidget(BoxLayout):
    def open(self, path, filename):
        with open(os.path.join(path, filename[0])) as f:
            print f.read()

    def selected(self, filename):
        print "selected: %s" % filename[0]


class Screen1App(Screen):
    myslider = ObjectProperty(None)
    user_request = ObjectProperty(None)
    request_terms = []

    def test_user_request(self):
        print self.user_request.text
        self.manager.current = 'Screen2'


    def run_alpha(self):
        """Alpha version of the run function
        use bibot script function.
        - A lot to do to parse nicely the request
        - create a label status
        """

        ## Parse the request
        request = self.user_request.text
        request = unidecode.unidecode(request) ## virer les accents
        request_terms = bibot.check_request_terms(request)

        if(request_terms == "NA"):
            ## TODO : display information in a label status
            print "Can't read input"
        else:

            self.request_terms = request_terms

            ## Run the process and go to load screen
            self.manager.current = 'Screen2'
            #time.sleep(2)
            #s2 = self.manager.get_screen('Screen2')
            #s2.display_someting(request_terms)
            
        



    ## Manage Files
	def open(self, path, filename):
		with open(os.path.join(path, filename[0])) as f:
			print f.read()

	def selected(self, filename):
		print "selected: %s" % filename[0]



class Screen2App(Screen):

    ## Screen 2 objects
    progress_bar = ObjectProperty(None)
    result_button = ObjectProperty(None)
    status_message = ObjectProperty(None)
    joke_message = ObjectProperty(None)
    run_button = ObjectProperty(None)


    def run_articles_selection(self):
        """
        The main function of bibot, fetch and select
        articles from pubmed based on their abstract.
        Functions within this function are called
        from the bibot script.
        """

        ## get request from Screen 1
        s1 = self.manager.get_screen('Screen1')
        request_term = s1.request_terms

        ## Dispaly Run information
        print "[INFO] PREPARE FOR RUN"
        self.status_message.text = "Prepare for run"

        ## Clean absract and meta folder
        print "[INFO] Cleaning directories"
        self.status_message.text = "Cleaning directories"
        for abstract_file in glob.glob("abstract/*.txt"):
            os.remove(abstract_file)
        for meta_data in glob.glob("meta/*.csv"):
            os.remove(meta_data)

        ## variables and file initialisation
        print "[INFO] Initialize log file"
        self.status_message.text = "Initialize log file"
        log_file = open("bibot.log", "w")

        ## Save request term in log file
        print "[INFO] TERMS : "
        request_line = ""
        for item in request_term:
            print "[-]\t"+str(item)
            request_line += str(item)+";"

        request_line = request_line[:-1]
        log_file.write(request_line+"\n")

        ## First interrogation of medline, get a
        ## huge list of articles possibly relevant and
        ## write the numbers of article in lig file
        print "[INFO] Large screening"
        big_list = bibot.get_huge_list_of_artciles(request_term)
        Total_number_of_articles = len(big_list)
        log_file.write("Total_number_of_articles;"+str(Total_number_of_articles)+"\n")
        print "[INFO] "+str(Total_number_of_articles) +" articles found"
        self.status_message.text = "Found "+str(Total_number_of_articles) +" articles to read"

        ## Test each articles retrieved from their pmid
        fetched = 0
        first_fiter_passed = 0
        last_filter_passed = 0
        cmpt = 0

        ## Start chrono
        long_request = False
        no_joke_display = True
        time_treshold = 20
        timer_a = datetime.datetime.now()


        for article in big_list:

            ## try to evaluate the article
            ## require a connection to the
            ## NCBI Server, if succed go on, 
            ## else wait 5 seconds and try again
            
            article_is_evaluated = False
            while(not article_is_evaluated):
                try:
                    valid = bibot.evaluate_article(article)
                    article_is_evaluated = True
                except:
                    print "[INFO] => CAN'T REACH NCBI, WAIT FOR 5 SECONDS "
                    self.status_message.text = "Looks like NCBI thinks you'are a hacker, I'm on it ;)"
                    now = datetime.datetime.now()
                    time_tag = str(now.hour)+"h:"+str(now.minute)+"m:"+str(now.day)+":"+str(now.month)
                    log_file.write("["+str(time_tag)+"];can't reach NCBI, wait for 5 seconds\n")
                    time.sleep(5)

            ## Get value for filters status
            filter_1_status = "FAILED"
            filter_2_status = "FAILED"
            if(valid[0]):
                fetched += 1
            if(valid[1]):
                first_fiter_passed += 1
                filter_1_status = "PASSED"
            if(valid[2]):
                last_filter_passed += 1
                filter_2_status = "PASSED"
            cmpt += 1

            ## Display progress on screen
            ## and status for each pmid in log file
            print "[RUN] => "+str(cmpt) +" [PROCESSED] || "+ str(fetched) + " [SELECTED] || FIRST FILTERS ["+filter_1_status+ "] || LAST FILTER ["+filter_2_status+ "] || "+str(float((float(cmpt)/float(Total_number_of_articles))*100)) + "% [COMPLETE]"
            log_file.write(">"+str(article)+";First_Filter="+str(filter_1_status)+";Last_filter="+str(filter_2_status)+"\n")

            ## update progress bar value
            progress = float((float(cmpt)/float(Total_number_of_articles))*100)
            self.progress_bar.value = progress
            display_progress = round(progress, 2)
            self.status_message.text = str(cmpt) +" articles processed, "+str(display_progress)+ "%"

            ## Check time spend on the task (i.e in the loop)
            ## if it take more than the time_treshold variable
            ## and the progress bar is under 50% the programm
            ## will start to display a few friendly message, 
            ## randomly selected in a list.
            timer_b = datetime.datetime.now()
            c = timer_b - timer_a
            c = divmod(c.days * 86400 + c.seconds, 60)
            if(c[1] > time_treshold and progress <= 50.0):
                long_request = True

            ## Select and display a joke
            if(long_request):
                joke_list = ["You should take a coffe", "In an other life, I'm so fast", "I bet I'm still reading faster than you !", "Boring ...", "You know that drinking tea is good for you ?", "I'm on it"]
                if(no_joke_display):
                    self.joke_message.text = joke_list[random.randint(0,len(joke_list)-1)]
                    no_joke_display = False
                else:
                    if(random.randint(0,100) > 80):
                        self.joke_message.text = joke_list[random.randint(0,len(joke_list)-1)]

        ## close log file
        log_file.close()

        ## Save the results files
        ## and folders
        now = datetime.datetime.now()
        time_tag = str(now.hour)+"h:"+str(now.minute)+"m:"+str(now.day)+":"+str(now.month)
        abstract_destination  = "SAVE/run_"+str(time_tag)+"/abstract"
        meta_destination = "SAVE/run_"+str(time_tag)+"/meta"
        log_destination = "SAVE/run_"+str(time_tag)+"/bibot.log"
        shutil.copytree("abstract", abstract_destination)
        shutil.copytree("meta", meta_destination)
        shutil.copy("bibot.log", log_destination)


        ## Generate results
        self.joke_message.text = "Computing results"
        bibot.plot_publications_years("meta")
        bibot.plot_country_stats("meta")


        ## Test
        self.run_button.size_hint_y = 0
        self.run_button.size_hint_x = 0
        self.run_button.text = ""
        self.run_button.opacity = 0
        self.run_button.disabled = True

        ## Go to result button
        self.result_button.size_hint_y = 2.5
        self.result_button.size_hint_x = 0.5
        self.result_button.text = "Go to Results"
        self.result_button.opacity = 1
        self.result_button.disabled = False


        ## Display end message
        self.joke_message.text = "Run completed"

            

    def run_articles_selection_wrapper(self):
        """
        Wrapper for the run articles_selection_function,
        design to call the target function in a Thread and
        therefore update the progress bar.
        """
        mythread = threading.Thread(target=self.run_articles_selection)
        mythread.start()









class Screen3App(Screen):
	pass


# Create the screen manager
sm = ScreenManager()
sm.add_widget(Screen1App(name='Screen1'))
sm.add_widget(Screen2App(name='Screen2'))
sm.add_widget(Screen3App(name='Screen3'))


class TestApp(App):
    title = "BIBOT"

    def build(self):

    	return sm
    	"""
        menu = Screen1App()
        #menu.test_user_request()
        print "value =",menu.myslider.value # <---- value here
        print "min =",menu.myslider.min # <--- min here
        return menu
		"""

if __name__ == "__main__":
    TestApp().run()
