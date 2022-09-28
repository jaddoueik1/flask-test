import json
from flask import Flask, render_template, request, url_for, redirect
from mypackage.economie_1 import economie
from mypackage.eco_affaire import eco_affaire
from mypackage.eco_affaire_er import eco_affaire_er
from mypackage.eco_premiere import eco_premiere
from mypackage.MTOW_dig import MTOWdig
from mypackage.itineraire_de_vol import distance_plot,company_plot,drapeau,direct

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/avion.html', methods=['GET', 'POST'])
def avion():
    if request.method == 'POST':
        npax = request.form.get('npax','300')
        raction = request.form.get('raction','3500')
        coefremp = request.form.get('coefremp','0.8')
        mach = request.form.get('mach','0.85')
        allongV = request.form.get('allongV','10')
        travel_class = request.form.get('options','classeEco')
        
        
        return redirect(url_for('resultats', npax=npax, raction=raction,coefremp=coefremp,mach=mach,allongV=allongV, travel_class=travel_class))

    return render_template('avion.html')

@app.route('/lesPropiete.html')
def les_propriete():
    return render_template('lesPropiete.html')

@app.route('/resultats.html?npax=<npax>&raction=<raction>&coefremp=<coefremp>&mach=<mach>&allongV=<allongV>&travel_class=<travel_class>')
def resultats(npax, raction,coefremp,mach,allongV, travel_class):

    
    if travel_class == 'classeEco':
        
        
        economie_object = economie(int(npax),int(raction),float(coefremp),float(mach),float(allongV))
        #economie_object.terminer()
        FN0, Masse_carb, MTOW, MC, ZFW, v_croi, MZFW = economie_object.terminer()
        prix = economie_object.TOC(FN0, Masse_carb, MTOW, MC, ZFW, v_croi, MZFW)
        #bg_remov()
        return render_template('resultats.html', prix=prix)
    
    if travel_class == 'classeEcoaff':
        
        mtow_object=MTOWdig(int(npax),int(raction),float(coefremp),float(mach),float(allongV))
        Sref,FN0, Masse_carb, MTOW, MC, ZFW, v_croi, MZFW=mtow_object.terminer()
       
        economie_aff_object = eco_affaire(int(npax),int(raction),float(coefremp),float(mach),float(allongV),float(Sref))
        economie_aff_object.terminer()
        prix = economie_aff_object.TOC(FN0, Masse_carb, MTOW, MC, ZFW, v_croi, MZFW)
        return render_template('resultats.html', prix=prix)
    
    if travel_class == 'classeEcofirstaff':
        
        mtow_object=MTOWdig(int(npax),int(raction),float(coefremp),float(mach),float(allongV))
        Sref,FN0, Masse_carb, MTOW, MC, ZFW, v_croi, MZFW=mtow_object.terminer()
        economie_aff_er_object = eco_affaire_er(int(npax),int(raction),float(coefremp),float(mach),float(allongV),float(Sref))
        economie_aff_er_object.terminer()
        prix = economie_aff_er_object.TOC(FN0, Masse_carb, MTOW, MC, ZFW, v_croi, MZFW)
        return render_template('resultats.html', prix=prix)
    
    if travel_class == 'classeEcofirst':
        
        mtow_object=MTOWdig(int(npax),int(raction),float(coefremp),float(mach),float(allongV))
        Sref,FN0, Masse_carb, MTOW, MC, ZFW, v_croi, MZFW=mtow_object.terminer()
        economie_er_object = eco_premiere(int(npax),int(raction),float(coefremp),float(mach),float(allongV),float(Sref))
        economie_er_object.terminer()
        prix = economie_er_object.TOC(FN0, Masse_carb, MTOW, MC, ZFW, v_croi, MZFW)
        return render_template('resultats.html', prix=prix)
        
    else:
        return render_template('avion.html')

@app.route('/route.html', methods=['GET', 'POST'])
def route():
    if request.method == 'POST':
        try:

            air_dep = request.form['depart']
            air_arrive = request.form['arrive']
            
        
            return redirect(url_for('trace', air_dep=air_dep, air_arrive=air_arrive))
            
        
        except:

            company=request.form['compagnie']
            return redirect(url_for('trace_company', company=company))



    else:
        
        return render_template('route.html')

@app.route('/resultatroute.html?air_dep=<air_dep>&air_arrive=<air_arrive>')
def trace(air_dep,air_arrive):

    distance=distance_plot(air_dep,air_arrive)
    pays_dep=drapeau(air_dep,air_arrive)[0]
    airport_dep=drapeau(air_dep,air_arrive)[1]
    pays_arr=drapeau(air_dep,air_arrive)[2]
    airport_arr=drapeau(air_dep,air_arrive)[3]
    vol=direct(air_dep,air_arrive)[0]
    nb_company=direct(air_dep,air_arrive)[1]

    return render_template('resultatroute.html', nb_company=nb_company,pays_dep=pays_dep,airport_dep=airport_dep, pays_arr=pays_arr,airport_arr=airport_arr,vol=vol,distance=distance,air_dep=air_dep,air_arrive=air_arrive)


@app.route('/resultat_route2.html?company=<company>')
def trace_company(company):

    company_plot(company)
    return render_template('resultat_route2.html')

    
if __name__ == '__main__':
    app.run(debug=False,host='0.0.0.0')

    