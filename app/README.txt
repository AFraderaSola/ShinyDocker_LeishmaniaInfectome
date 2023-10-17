##################### Leish infectome shiny app
##### Scripts

# myui.R
User interface script. Here you deine how the web app will look and what are the user's inputs.

# myuiTabs.R
User interface script. Here you deine how the web app will look and what are the user's inputs. This one looks nicer, with tabs.

# myserver.R
Functions that shape the data, plot the results and so on. 

# app.R
Script for running the shiny app. It just sources the myui and myserver scripts and generates the web app.

##### InputFiles

# QuantifiedProteins_.*.csv
Quantified proteins for each leishmania species after LFQ analysis

# PrepareDataShiny.R
It generates the database and test_data dataframes (in .csv) needed for running the app
