############################################# Dossier Optimisation - Recherche Operationnelle ###################### 
#PERENON Clement, REGAZZETTI Lea
#M1 Informatique
#Promotion 2020-2021

############################################### ALGORITHME DE PRIM #################################################
Prim=function(X,A)
{
  #Fonction qui determine un arbre partiel de poids minimal
  #INPUT:
  #X l'ensemble des sommets
  #A est la matrice d'adjacence ponderee
  #OUTPUT:
  #U_prim la liste des aretes de l'arbre partiel de poids minimal
  i=1 #on choisit au hasard un sommet de X
  x_prim=c(i) #ce sommet choisi est integre dans la liste des sommets de l'arbre partiel de poids minimal
  u_prim=c() #initialisation de la liste des aretes de l'arbre partiel de poids minimal
  cpt=1 #compteur des iterations
  
  while(setequal(X,x_prim)==FALSE) #tant que tous les sommets X du graphe ne sont pas contenus dans les sommets X' de l'arbre partiel de poids minimal
  {
    j=intersect(X,x_prim) #Sommet j appartennat a la liste des sommets de l'arbre partiel de poids minimal
    k=setdiff(X,x_prim) #Sommet k n'appartennat pas a la liste des sommets de l'arbre partiel de poids minimal
    
    if (length(j)==1 || length(k)==1) 
    {
      #si il n'y a qu'un sommet dans la liste des sommets, on transformera sous forme de matrice
      successeurs=which(matrix(A[j,k]!=0,nrow=length(j),ncol=length(k)),arr.ind=TRUE) #indice des successeurs de j parmi les sommets k 
    }
    else
    {
      successeurs=which(A[j,k]!=0,arr.ind=TRUE) #indice des successeurs de j parmi les sommets k
    }
    
    successeurs=cbind(j[successeurs[,1]],k[successeurs[,2]]) #indices des successeurs de j parmi les sommets k dans la matrice d'adjacence ponderee A
    indice=which.min(A[successeurs]) #indice du successeur ayant le poids minimal dans la matrice "successeurs"
    res=cbind(successeurs[indice,1],successeurs[indice,2]) #indices du successeur ayant le poids minimal dans la matrice A
    j=res[1] #sommet appartenant a la liste des sommets de l'arbre partiel ayant permis de determiner le sommet k
    k=res[2] #sommet ayant le poids minimal n'appartenant pas a la liste des sommets de l'arbre partiel
    arete=c(j,k) #indice des sommets j et k constituant l'arete de poids minimal
    x_prim=union(x_prim,k) #liste des sommets de l'arbre partiel de poids minimal
    u_prim[[cpt]]=arete #liste des aretes de l'arbre partiel de poids minimal
    cpt=cpt+1
  }
  
  print('X prim : ')
  print(x_prim)
  
  print('U prim :')
  print(u_prim)
}


#Exemple du projet : 
X=c(1:7) 
A=rbind(c(0,5,8,0,0,0,0),c(5,0,0,4,2,0,0),c(8,0,0,0,5,2,0),c(0,4,0,0,0,0,7),c(0,2,5,0,0,0,3),c(0,0,2,0,0,0,3),c(0,0,0,7,3,3,0))
print(A)
Prim(X,A)

#Exemple du cours: 
X=c(1:5)
A=rbind(c(0,6,5,1,10),c(6,0,4,2,8),c(5,4,0,7,0),c(1,2,7,0,0),c(10,8,0,0,0))
print(A)
Prim(X,A)

X=c(1:6)
A=rbind(c(0,1,4,2,0,0),c(1,0,8,5,5,0),c(4,8,0,4,0,6),c(2,3,4,0,0,5),c(0,5,0,0,0,2),c(0,0,6,5,2,0))
print(A)
Prim(X,A)

############################################### ALGORITHME DE FORD BELLMAN ##########################################

Ford_Bellman = function(X,A,s)
{
  #Determine les longueurs des plus courts chemins entre un sommet s et les autres sommets
  #INPUT:
  #X l'ensemble des sommets
  #A est la matrice d'adjacence
  #s un sommet de X (sommet de depart)
  #OUTPUT:
  #N les longueurs des plus courts chemins entre s et tous les autres sommets de X

  p=c()
  p[s]=0 #Pour le sommet de depart, la longueur vaut 0 
  s_barre=setdiff(X,s) 
  for (i in s_barre)
  {
    p[i] = Inf #initialise la longueur de tous les autres sommets a l'infini
  }
  temp = p #variable temporaire contenant la valeur courante des longueurs 
  change = 1 #variable permettant de tester si une des valeurs de p[i] change
  
  print("Initialisation : ")
  print(p)
  
  cpt=1
  while (change > 0) #si la valeur de la variable "change" est superieure a O, cela signifie qu'une des valeurs de p[i] a change lors de l'iteration
  {
    print(paste("Iteration n°" , cpt))
    for (i in s_barre) #pour chaque sommet i different du sommet de depart
    {
      pre=which(A[,i]!=0) #indices des predecesseurs du sommet i
      res=min(p[pre]+A[pre,i]) #pour tous les predecesseurs de i on prend la longueur du chemin de s jusqu'au sommet predecesseur du sommet i + la longueur du sommet predecesseur jusqu'au sommet i et on recupere la longueur minimale
      p[i]=min(p[i],res) #valeur minimale entre le resultat obtenu precedemment et la valeur courante de p[i]
    }
    change=sum(as.numeric(temp!=p)) #somme des valeurs differentes entre temp et p est superieure a 0
    temp = p #mise a jour de la variable temporaire
    cpt=cpt+1
    print(p) #affichage de la valeur courante du vecteur p
  }
}

#Exemple du cours
A=rbind(c(0,7,8,0,0,0),c(0,0,0,4,1,2),c(0,2,0,0,0,2),c(0,0,0,0,0,0),c(0,0,-2,2,0,0),c(0,0,0,0,3,0))
A
X=1:6
s=1
Ford_Bellman(X,A,s)

#Exemple 1 du projet: 
X=c(1:7) 
A=rbind(c(0,5,8,0,0,0,0),c(5,0,0,4,-2,0,0),c(8,0,0,0,5,2,0),c(0,4,0,0,0,0,7),c(0,-2,5,0,0,0,3),c(0,0,2,0,0,0,-3),c(0,0,0,7,3,-3,0))
print(A)
s=1
Ford_Bellman(X,A,s)

#Exemple 2 du projet: 
X=c(1:7) 
A=rbind(c(0,5,8,0,0,0,0),c(5,0,0,4,0,0,0),c(8,0,0,0,5,2,0),c(0,4,0,0,0,0,0),c(0,0,5,0,0,0,3),c(0,0,2,0,0,0,-3),c(0,0,0,0,3,-3,0))
print(A)
s=7
Ford_Bellman(X,A,s)




############################################ ALGORITHME DE FORD FULKERSON ###########################################

Ford_Fulkerson = function(X,A,s,p)
{
  n = nrow(A) #Nombre de ligne de la matrice A en parametre
  sommets = X #ensemble de sommet
  flot_debut = matrix(0,nrow=n,ncol=n) 
  marque = matrix(0,ncol=3,nrow=n) #Je définit une matrice pour les 3 paramètres de ms
  s_marque = s #affecte le sommet source aux sommets marqués (car pas besoin de le modif)
  s_pasmarque = setdiff(sommets,s_marque) #On met les sommets qui pas marqués pour le moment
  matrice_arete_bool = (A-flot_debut>0) #matrice de booléen des arêtes
  indice_aretes = which(matrix(matrice_arete_bool[s_marque,s_pasmarque]==TRUE,nrow=length(s_marque),ncol=length(s_pasmarque)), arr.ind=TRUE) #le which nous permet de récupérer les indices des aretes
  flot_maximum = 0 #Initialiser le flot maximum
  
  ##########################################################################
  
  while (TRUE)
  {
    #On cherche si deux sommets satisfont les conditions
    if(length(indice_aretes>0))
    {
      i=s_marque[indice_aretes[1,1]]
      j=s_pasmarque[indice_aretes[1,2]]
    }
    else
    {
      break;
    }
    #Si Cij - Phiij > 0 (Cij capacité de l'arc et Phiij le flot) alors :
    if (A[i,j]-flot_debut[i,j] > 0)
    {
      #Calculer la valeur maximum et on marque avec 1 qui correspondra au troisième argument, +
      a = abs(min(marque[i,2],A[i,j]-flot_debut[i,j])) 
      marque[j,] = c(i,a,1)
    }
    else
    {
      if (flot_debut[j,i] >-1)
      { 
        #Marquer avec -1 qui correspondera au troisième argument, -
        a = abs(min(marque[i,2],flot_debut[j,i]))
        marque[j,] = c(i,a,-1)
      }
    }
    #Réinitialiser la matrice des indice_aretes
    s_marque = union(s_marque,j)#On fait l'union des sommets marqués et j 
    s_pasmarque = setdiff(s_pasmarque,s_marque)#On récupère ceux qui ne sont pas marqués ici
    indice_aretes = which(matrix(matrice_arete_bool[s_marque,s_pasmarque]==TRUE,nrow=length(s_marque),ncol=length(s_pasmarque)),arr.ind=TRUE)
    
    ######################################################################
    
    if(j==p) #Si notre j correspond au puit, alors :
    {
      #On calcule l'accroissement du flot initla
      if(is.element(p,s_marque))
      { 
        while(j!=s) #Si j est différent de la source alors :
        {
          if(marque[j,3] == 1) #Si mj est +, alors :
          {
            #Modifier le flot entre i et j
            flot_debut[marque[j,1],j] = flot_debut[marque[j,1],j]+marque[p,2] 
          }
          else if(marque[j,3] == -1)
          {
            #Modifier le flot entre j et i
            flot_debut[j,marque[j,1]] = flot_debut[j,marque[j,1]]-marque[p,2] 
          }
          j = marque[j,1] 
        }
        #Nouvelle valeur du flot maximum
        flot_maximum = flot_maximum+marque[p,2]
      }
      marque[s,c(2,3)] = c(Inf,1) 
      s_marque = s 
      s_pasmarque = setdiff(sommets,s_marque)
      matrice_arete_bool = (A-flot_debut>0)
      indice_aretes = which(matrix(matrice_arete_bool[s_marque,s_pasmarque]==TRUE,nrow=length(s_marque),ncol=length(s_pasmarque)),arr.ind=TRUE)
    }
  }
  
  #########################################################################
  
  #Calculer l'accroissement de la valeur du flot réalisable
  if(is.element(p,s_marque))
  { 
    while(j != s)
    {
      if(marque[j,3] == 1)
      {
        #Modifier le flot entre i et j
        flot_debut[marque[j,1],j] = flot_debut[marque[j,1],j]+marque[p,2]  
      }
      else if(marque[j,3] == -1)
      {
        #Modifier le flot entre j et i
        flot_debut[j,marque[j,1]] = flot_debut[j,marque[j,1]]-marque[p,2] 
      }
      j = marque[j,1] 
    }
  }
  
  #######################################################################
  else
  { 
    return(flot_maximum) #on retourne le flot maximum
  }
  
}

#Exemple du projet
A = rbind(c(0,1,0,0,1,0),c(0,0,0,0,0,0),c(0,0,0,0,0,0),c(0,4,4,0,0,4),c(0,5,5,0,0,0),c(6,0,6,0,0,0))
X=1:nrow(A)
Ford_Fulkerson(X,A,4,2)

#Exemple du cours
A=rbind(c(0,5,8,0,0,0,0),c(0,0,0,4,2,0,0),c(0,0,0,0,5,2,0),c(0,0,0,0,0,0,7),c(0,0,0,0,0,0,3),c(0,0,0,0,0,0,3),c(0,0,0,0,0,0,0))
X=1:nrow(A)
Ford_Fulkerson(X,A,1,7)

#Partiel 2017
#flot_debut = rbind(c(0,5,5,0,0,0,0),c(0,0,0,0,5,0,0),c(0,0,0,5,0,0,0),c(0,0,0,0,0,0,0),c(0,0,0,0,0,0,10),c(0,0,0,0,0,0,0),c(0,0,0,0,0,0,0))
A=rbind(c(0,6,8,0,0,0,0),c(0,0,8,0,5,0,0),c(0,0,0,6,10,2,0),c(0,0,0,0,0,5,0),c(0,0,0,0,0,0,10),c(0,0,2,0,0,0,4),c(0,0,0,0,0,0,0))
X=1:nrow(A)
Ford_Fulkerson(X,A,1,7)

#Partiel 2020
A=rbind(c(0,10,2,5,0,0,0),c(0,0,5,0,3,0,0),c(0,0,0,4,0,0,0),c(0,0,0,0,2,10,0),c(0,0,0,0,0,0,7),c(0,0,0,0,0,0,8),c(0,0,0,0,0,0,0))
X=1:nrow(A)
Ford_Fulkerson(X,A,1,7)