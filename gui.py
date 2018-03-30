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
import datetime

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
            text: "Screen 2"
            color: 0, 0, 0, 1
            outline_color: 0, 0.5, 0.5, 1
            font_size: 30
        ProgressBar:
        	id: pb
        	size_hint_x: .5
        	size_hint_y: None
        	height: '48dp'
			value: 0
		Button
			size_hint: 0.3, 0.4
			text: "Run"
			on_release: root.play_with_pb_values()
            on_release: root.display_something()
		Button
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
	            text: "Screen 3 - Results 1"
	            color: 0, 0, 0, 1
	            outline_color: 0, 0.5, 0.5, 1
	            font_size: 30

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
    progress_bar = ObjectProperty(None)
    result_button = ObjectProperty(None)



    def play_with_pb_values(self):
        print self.progress_bar.value

        if(self.progress_bar.value < 100):
			self.progress_bar.value += 25
        else:
            print "Display Results"
            self.result_button.size_hint_y = .5
            self.result_button.size_hint_x = .5
            self.result_button.text = "Go to Results"
            self.result_button.opacity = 1
            self.result_button.disabled = False

    
    def display_something(self):
        """Test an adapatation of run function
        from BIBOT
        """

        ## get request from Screen 1
        s1 = self.manager.get_screen('Screen1')
        request_term = s1.request_terms

        ## Dispaly Run information
        print "[INFO] PREPARE FOR RUN"


        ## Clean absract and meta folder
        print "[INFO] Cleaning directories"
        for abstract_file in glob.glob("abstract/*.txt"):
            os.remove(abstract_file)
        for meta_data in glob.glob("meta/*.csv"):
            os.remove(meta_data)

        ## variables and file initialisation
        print "[INFO] Initialize log file"
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

        ## Test each articles retrieved from their pmid
        fetched = 0
        first_fiter_passed = 0
        last_filter_passed = 0
        cmpt = 0

        for article in big_list:

            ## try to evaluate the article
            ## require a connection to the
            ## NCBI Server, if succed go on, 
            ## else wait 5 seconds and try again

            article_is_evaluated = False
            while(not article_is_evaluated):
                try:
                    #print "|| TRY TO PROCESS ARTICLE "+str(article)+ " ||"
                    valid = bibot.evaluate_article(article)
                    article_is_evaluated = True
                except:
                    #print "|| CAN'T REACH NCBI, WAIT FOR 5 SECONDS ||"
                    print "[INFO] => CAN'T REACH NCBI, WAIT FOR 5 SECONDS "
                    now = datetime.datetime.now()
                    time_tag = str(now.hour)+"h:"+str(now.minute)+"m:"+str(now.day)+":"+str(now.month)
                    log_file.write("["+str(time_tag)+"];can't reach NCBI, wait for 5 seconds\n")
                    time.sleep(5)


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

            ## update progress bar value -> not working for now
            self.progress_bar.value += 10

            print "[PROGRESS] => "+str(self.progress_bar.value)

            #float((float(cmpt)/float(Total_number_of_articles))*100)
            


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
