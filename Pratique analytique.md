**Chargement des packages nécessaires pour le scipt.**

```{r}
library(RMark)
library(tidyverse)
```

# Sterne Dougall (Sterna dougallii), les plus grosses colonies européennes : Royaume Uni et Irelande.

## I - Importation et formatage des données pour le script :

```{r}
setwd("C:/Users/UTILISATEUR/Documents/cours/M1/S8/Analyses de séquences/Modélo cam/rapport")
Stern <- read.table('ROSEATE_METAPOP.TXT')
Stern


names(Stern) = "ch"
Stern$freq = 1
Stern$ch = as.character(Stern$ch)
head(Stern)
```

## II - Création du "process data" et du "design data" :

Dans le process data on inclue l'année de début des mesures et le site originaire des sternes.

```{r}
dp = process.data(Stern, model = "Multistrata", strata.labels = c("C","L", "R"), begin.time = 1992)
dp
```

On construit ensuite le design data à partir du process data créé ultérieurement.

```{r}
ddl = make.design.data(dp) 
ddl
```

## III - Fabrication du premier modèle pour p (probabilité d'observation) :

Dans un premier temps on essaie de trouver le meilleur modèle pour p comme il nous a été conseillé. Par la suite, on modifie le ddl pour y inclure des classes d'âges comme recommandé par les auteurs. Mais aussi, les distances entre les différents sites pour le paramètre psi.

P permet de connaitre la probabilité d'observation ainsi on peut savoir si les oiseaux ont été observés régulièrement et ainsi connaitre la significativité des observations.

Ainsi, pour trouver le meilleur modèle on essaie d'inclure différents paramètres tel que le "time" ou le "ageclass" etc...

```{r}
run.Stern=function()
{
  # Process data (dp)
  dp = process.data(Stern, model = "Multistrata", strata.labels = c("C","L", "R"), begin.time = 1992)
  # Design data (ddl)
  ddl = add.design.data(dp, ddl, parameter = "S", type = "age", bins=c(0,2,25),name="ageclass")
  ddl = add.design.data(dp, ddl, parameter = "Psi", type = "age", bins=c(0,2,25),name="ageclass",replace = TRUE)
  ddl = add.design.data(dp, ddl, parameter = "p", type = "age", bins=c(0,2,25),name="ageclass",replace = TRUE)
  ddl$Psi$distance=0
  ddl$Psi$distances[ddl$Psi$stratum=="C"&ddl$Psi$tostratum=="L"]=525
  ddl$Psi$distances[ddl$Psi$stratum=="C"&ddl$Psi$tostratum=="R"]=350
  ddl$Psi$distances[ddl$Psi$stratum=="L"&ddl$Psi$tostratum=="R"]=200
  ddl$Psi$distances[ddl$Psi$stratum=="L"&ddl$Psi$tostratum=="C"]=525
  ddl$Psi$distances[ddl$Psi$stratum=="R"&ddl$Psi$tostratum=="C"]=350
  ddl$Psi$distances[ddl$Psi$stratum=="R"&ddl$Psi$tostratum=="L"]=200
  # Formules pour trouver le meilleur modèle pour p
  Psi.distance=list(formula=~distance)
  p.stratum=list(formula=~stratum)
  p.stratum.time=list(formula=~stratum+time)
  p.stratum.ac.time=list(formula=~stratum+time+ageclass)
  p.stratumtimeac=list(formula=~stratum*time*ageclass)
  p.stratumac=list(formula=~stratum*ageclass)
  p.stratumtime=list(formula=~stratum*time)
  p.stratum.Time=list(formula=~stratum+Time)
  p.stratum.ac=list(formula=~stratum+ageclass)
  p.stratum.ac.Time=list(formula=~stratum+Time+ageclass)
  p.dot=list(formula=~1)
  S.stratum=list(formula=~stratum)
  model.list=create.model.list("Multistrata")
  Stern.results=mark.wrapper(model.list,data=dp,ddl=ddl)
  return(Stern.results)
}
Stern.results=run.Stern()
Stern.results
```

Le meilleur modèle obtenu pour p est donc "p.stratumtimeac". Il comprend donc le temps continu, la classe d'age et le site d'observation. On sait que c'est le meilleur modèle car c'est le modèle avec le plus faible AICc (100011.4 ) et son poids de 1.

## IV - Fabrication du second modèle pour Phi (survie) :

Le second modèle a être recherché est le modèle pour la survie. Pour ce modèle comme le précédent on conservera le même dp et ddl. De plus, on conservera aussi le meilleur modèle p pour continuer à travailler.

Encore une fois, pour ce modèle ce qui est recherché est l'AICc le plus faible ainsi on va modifié différents paramètres du modèle de survie S pour obtenir le meilleur.

```{r}
run.Stern1=function()
{
  # Process data (dp)
  dp = process.data(Stern, model = "Multistrata", strata.labels = c("C","L", "R"), begin.time = 1992)
  # Design data (ddl)
  ddl = add.design.data(dp, ddl, parameter = "S", type = "age", bins=c(0,2,25),name="ageclass")
  ddl = add.design.data(dp, ddl, parameter = "Psi", type = "age", bins=c(0,2,25),name="ageclass",replace = TRUE)
  ddl = add.design.data(dp, ddl, parameter = "p", type = "age", bins=c(0,2,25),name="ageclass",replace = TRUE)
  ddl$Psi$distance=0
  ddl$Psi$distances[ddl$Psi$stratum=="C"&ddl$Psi$tostratum=="L"]=525
  ddl$Psi$distances[ddl$Psi$stratum=="C"&ddl$Psi$tostratum=="R"]=350
  ddl$Psi$distances[ddl$Psi$stratum=="L"&ddl$Psi$tostratum=="R"]=200
  ddl$Psi$distances[ddl$Psi$stratum=="L"&ddl$Psi$tostratum=="C"]=525
  ddl$Psi$distances[ddl$Psi$stratum=="R"&ddl$Psi$tostratum=="C"]=350
  ddl$Psi$distances[ddl$Psi$stratum=="R"&ddl$Psi$tostratum=="L"]=200
  # Formules pour trouver le meilleur le modèle de survie
  Psi.distance=list(formula=~distance)
  p.stratum.ac.time=list(formula=~stratum*time*ageclass)
  S.stratum=list(formula=~stratum)
  S.stratum.ac=list(formula=~stratum+ageclass)
  S.stratum.time=list(formula=~stratum+time)
  S.stratum.time.ac=list(formula=~stratum+time+ageclass)
  S.stratum.Time=list(formula=~stratum+Time)
  S.stratum.Time.ac=list(formula=~stratum+Time+ageclass)
  S.stratumtimeac=list(formula=~stratum*time*ageclass)
  S.stratumtime=list(formula=~stratum*time)
  S.stratumac=list(formula=~stratum*ageclass)
  model.list=create.model.list("Multistrata")
  Stern.results1=mark.wrapper(model.list,data=dp,ddl=ddl)
  return(Stern.results1)
}
Stern.results1=run.Stern1()
Stern.results1
```

Comme pour le modèle précéden,t on retrouve le meilleur modèle d'AICc avec les paramètres de temps continu, la classe d'age et le site d'observation. C'est donc le modèle "S.stratumtimeac" ce qui permet de faire baisser l'AICc à 97360.04 et se sera le seul conservé car il possède un weight de 1.

## V - Fabrication du troisième modèle pour Psi (probabilité de passer d'une colonie à une autre) :

Ce dernier modèle conserne Psi c'est à dire la probabilité qu'une sternie passe d'une colonie à une autre. Pour ce modèle comme pour les précédents on conservera le même dp et ddl. De plus, on conservera aussi le meilleur modèle p et S pour continuer à travailler.

Encore une fois, pour ce modèle ce qui est recherché est l'AICc le plus faible ainsi on va modifier différents paramètres du modèle Psi pour essayer d'obtenir le meilleur possible.

Cette fois si on a décidé par sécurité de diviser la recherche en deux fonctions. Celles-ci étant très lourde pour l'ordinateur nous avons décidé de choisir la sécurité. Mais cela n'empèche en rien la recherche du meilleur modèle entre chaques fonctions.

Donc cette première fonction nous avons décidé de travailler sur le paramètre distance que nous avons rajouté dans le ddl.

```{r}
run.Stern2=function()
{
  # Process data (dp)
  dp = process.data(Stern, model = "Multistrata", strata.labels = c("C","L", "R"), begin.time = 1992)
  # Design data (ddl)
  ddl = add.design.data(dp, ddl, parameter = "S", type = "age", bins=c(0,2,25),name="ageclass")
  ddl = add.design.data(dp, ddl, parameter = "Psi", type = "age", bins=c(0,2,25),name="ageclass",replace = TRUE)
  ddl = add.design.data(dp, ddl, parameter = "p", type = "age", bins=c(0,2,25),name="ageclass",replace = TRUE)
  ddl$Psi$distance=0
  ddl$Psi$distances[ddl$Psi$stratum=="C"&ddl$Psi$tostratum=="L"]=525
  ddl$Psi$distances[ddl$Psi$stratum=="C"&ddl$Psi$tostratum=="R"]=350
  ddl$Psi$distances[ddl$Psi$stratum=="L"&ddl$Psi$tostratum=="R"]=200
  ddl$Psi$distances[ddl$Psi$stratum=="L"&ddl$Psi$tostratum=="C"]=525
  ddl$Psi$distances[ddl$Psi$stratum=="R"&ddl$Psi$tostratum=="C"]=350
  ddl$Psi$distances[ddl$Psi$stratum=="R"&ddl$Psi$tostratum=="L"]=200
  # Formules pour trouver le meilleur modèle pour psi
  Psi.distance=list(formula=~distance)
  Psi.distance.time=list(formula=~distance+time)
  Psi.distance.Time=list(formula=~-distance+Time)
  Psi.distance.ac=list(formula=~distance+ageclass)
  Psi.distance.time.ac=list(formula=~distance+time+ageclass)
  Psi.distance.Time.ac=list(formula=~distance+Time+ageclass)
  Psi.r=list(formula=~ -1+stratum:tostratum)
  Psi.distance.co=list(formula=~distance+cohort)
  p.stratum.ac.time=list(formula=~stratum*time*ageclass)
  S.stratum.time.ac=list(formula=~stratum*time*ageclass)
  model.list=create.model.list("Multistrata")
  Stern.results2=mark.wrapper(model.list,data=dp,ddl=ddl)
  return(Stern.results2)
}
Stern.results2=run.Stern2()
Stern.results2
```

La suite des recherches va se faire avec les paramètres stratum et tostratum.

```{r}
run.Stern3=function()
{
  # Process data (dp)
  dp = process.data(Stern, model = "Multistrata", strata.labels = c("C","L", "R"), begin.time = 1992)
  # Design data (ddl)
  ddl = add.design.data(dp, ddl, parameter = "S", type = "age", bins=c(0,2,25),name="ageclass")
  ddl = add.design.data(dp, ddl, parameter = "Psi", type = "age", bins=c(0,2,25),name="ageclass",replace = TRUE)
  ddl = add.design.data(dp, ddl, parameter = "p", type = "age", bins=c(0,2,25),name="ageclass",replace = TRUE)
  ddl$Psi$distance=0
  ddl$Psi$distances[ddl$Psi$stratum=="C"&ddl$Psi$tostratum=="L"]=525
  ddl$Psi$distances[ddl$Psi$stratum=="C"&ddl$Psi$tostratum=="R"]=350
  ddl$Psi$distances[ddl$Psi$stratum=="L"&ddl$Psi$tostratum=="R"]=200
  ddl$Psi$distances[ddl$Psi$stratum=="L"&ddl$Psi$tostratum=="C"]=525
  ddl$Psi$distances[ddl$Psi$stratum=="R"&ddl$Psi$tostratum=="C"]=350
  ddl$Psi$distances[ddl$Psi$stratum=="R"&ddl$Psi$tostratum=="L"]=200
  # Formules pour trouver le meilleur modèle pour psi
  Psi.r=list(formula=~ -1+stratum:tostratum)
  Psi.r.time=list(formula=~ -1+stratum:tostratum:time)
  Psi.r.time.ac=list(formula=~ -1+stratum:tostratum:time:ageclass)
  Psi.rplustimeplusac=list(formula=~ stratum+tostratum+time+ageclass)
  Psi.rplustimeplusacplusdistance=list(formula=~ stratum+tostratum+time+ageclass+distance)
  Psi.rplusdistanceplusdistance.time.ac=list(formula=~ -1+stratum+distance:tostratum+distance:time:ageclass)
  Psi.rplusdistance.time.ac=list(formula=~ -1+stratum+distance:tostratum:time:ageclass)
  Psi.rplustime=list(formula=~ stratum+tostratum+time)
  Psi.rplus=list(formula=~ stratum+tostratum)
  Psi.rr=list(formula=~ stratum*tostratum)
  Psi.rtime=list(formula=~ stratum*tostratum*time)
  Psi.rtimeac=list(formula=~ stratum*tostratum*time*ageclass)
  Psi.rtimeacdistance=list(formula=~ stratum*tostratum*time*ageclass*distance)
  p.stratum.ac.time=list(formula=~stratum*time*ageclass)
  S.stratum.time.ac=list(formula=~stratum*time*ageclass)
  model.list=create.model.list("Multistrata")
  Stern.results3=mark.wrapper(model.list,data=dp,ddl=ddl)
  return(Stern.results3)
}
Stern.results3=run.Stern3()
Stern.results3
Stern.results2
```

On peut maintenant comparer les différents modèles pour choisir le meilleur. L'AICc le plus faible obtenu est de 95708.64 avec un weight de 1.

Voici donc notre modèle finale comportant tous les meilleurs modèles trouvés auparavant :

Psi.r.time.ac=list(formula=\~ -1+stratum:tostratum:time:ageclass) p.stratumtimeac=list(formula=\~stratum\*time\*ageclass) S.stratumtimeac=list(formula=\~stratum\*time\*ageclass)

Il parait assez logique que se soit les modèles avec \* qui ont le résultat avec le plus faible AICc car "time" et "ageclass" sont en interaction avec le paramètre stratum.

```{r}
run.Sternf=function()
{
  # Process data (dp)
  dp = process.data(Stern, model = "Multistrata", strata.labels = c("C","L", "R"), begin.time = 1992)
  # Design data (ddl)
  ddl = add.design.data(dp, ddl, parameter = "S", type = "age", bins=c(0,2,25),name="ageclass")
  ddl = add.design.data(dp, ddl, parameter = "Psi", type = "age", bins=c(0,2,25),name="ageclass",replace = TRUE)
  ddl = add.design.data(dp, ddl, parameter = "p", type = "age", bins=c(0,2,25),name="ageclass",replace = TRUE)
  ddl$Psi$distance=0
  ddl$Psi$distances[ddl$Psi$stratum=="C"&ddl$Psi$tostratum=="L"]=525
  ddl$Psi$distances[ddl$Psi$stratum=="C"&ddl$Psi$tostratum=="R"]=350
  ddl$Psi$distances[ddl$Psi$stratum=="L"&ddl$Psi$tostratum=="R"]=200
  ddl$Psi$distances[ddl$Psi$stratum=="L"&ddl$Psi$tostratum=="C"]=525
  ddl$Psi$distances[ddl$Psi$stratum=="R"&ddl$Psi$tostratum=="C"]=350
  ddl$Psi$distances[ddl$Psi$stratum=="R"&ddl$Psi$tostratum=="L"]=200
  # Formule pour le modèle
  Psi.r.time.ac=list(formula=~ -1+stratum:tostratum:time:ageclass)
  p.stratum.ac.time=list(formula=~stratum*time*ageclass)
  S.stratum.time.ac=list(formula=~stratum*time*ageclass)
  model.list=create.model.list("Multistrata")
  Stern.resultsf=mark.wrapper(model.list,data=dp,ddl=ddl)
  return(Stern.resultsf)
}
Stern.resultsf=run.Sternf()
Stern.resultsf
summary(run.Sternf())


```

En revanche nous n'avons pas réussi à faire fonctionner le "model averaging" pour estimer les paramètres une erreur nous a bloqué lors de cette dernière étape.

L'AIC est un outil utilisé en statitistique pour comparer différents modèles et évaluer leur ajustement aux données. Ce modèle prouve bien que psi et phy, le taux de substitution et le taux de survie, impactent p, le taux de détection de la stern. De plus, le weight de 1 semble assez logique car l'AICc étant très élevé il y a donc plus de chance que le l'écart entre les modèles soient très significatif.

Ce modèle nous permet d'observer un taux de survie qui varie entre les années mais aussi entre les sites, de plus la classe d'âge influe aussi sur la survie. La survie est plus faible chez la classe d'age 1 (entre 0 et 2 ans) que sur la seconde classe d'âge, on observe que cette variation peut être de 50% entre les sternes juvéniles et les sternes aldutes. L'année semble aussi provoquer des variations sur la survie allant jusqu'à un écart de 30% chez les adultes. La localité aura aussi un impact car sur la même année et sur la même classe d'âge la probabilité de survie ne sera pas la même entre L, R ou encore C. En revanche, on observe que le paramètre distance ne permet pas d'obtenir un meilleur modèle et même au contraire cela le détériore en augmentant l'AICc final.

Le taux de survie va être impacté par le taux de subsitution car la survie varie d'une île à une autre de plus il semble y avoir un plus gros mouvement des sternes vers le site R et un très faible mouvement des sterns entre le site C et L.

Le model average n'ayant pas fonctionné il est compliqué d'obtenir d'information et de donc de donner plus de conclusion.
