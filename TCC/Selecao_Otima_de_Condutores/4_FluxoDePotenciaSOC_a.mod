#------------------------------#
# Fluxo de Potencia Nao Linear #
#------------------------------#

#-- Definir os Conjuntos --#

set Ob; #Conjunto das barras 
set Ol within Ob cross Ob; #Conjunto das linhas, que depende de Ob
set Oc;

#-- Definir os Parametros --#

# Dados das barras
param Tb{Ob};	#Tipo de barra (1: SE, 0: Carga)
param PD{Ob};	#Potencia ativa demandada
param QD{Ob};	#Potencia Reativa Demandada
param QC{Ob};	#Pot Reativa dos capacitores na barra

param Vnom;
param Vmin;
param Vmax;

# Dados das linhas
param Rc{Oc};	#Resistencia
param Xc{Oc};	#Reatancia indutiva
param Icmax{Oc};	#Corrente maxima
param Z2{Oc};	#Z^2 = R^2 + X^2
param Cc{Oc};
param L{Ol};
param Linea{Ol};
param TC{Ol};

#Dados condutores

param Kp;
param Ke;
param T;
param LSF;
param IDF;

# linearização

var Fi{Ol,Oc} >= 0;
param M1 = 10^14;
param M2 = 10^7;
param M3 = 10^7;
var Beta{Ol,Oc} >= 0;
var Gama{Ol,Oc} >= 0;



#-- Definir as Variaveis --#

# Vars das barrasS
var PS{Ob};		#Potencia ativa
var QS{Ob};		#Potencia reativa

# Vars das linhas
var P{Ol};		#Potencia ativa na linha
var Q{Ol};		#Potencia reativa na linha

# Tensoes e correntes
var Isqr{Ol} >= 0;		#corrente nas linhas
var Vsqr{Ob};		#tensões nas barras

#----------------------------------------

# Para a linearização Vsqr * Isqr

param S = 1;
param DeltaV = (Vmax^2-Vmin^2)/(S+1);
var xv{Ob, s in 1..S}, binary;
var Pc{Ob, s in 1..S};
var W{Ol,Oc}, binary;
var A{Ol, Oc, s in 1..S} >= 0,  binary;

#----------------------------------------

# Dados para a Linearizacao P^2 Q^2

param Y = 50;
param DS{Ol};
param ms{Ol, y in 1..Y};

var DP{Ol, y in 1..Y} >= 0;
var DQ{Ol, y in 1..Y} >= 0;

var Pmax{Ol}>= 0;
var Pmin{Ol}>= 0;
var Qmax{Ol}>= 0;
var Qmin{Ol}>= 0;
#----------------------------------------

#Variaveis para geração distribuda




#-- Funcao Objetivo --#

minimize FuncaoObjetivo: (Kp + (Ke + T + LSF)) * sum{(i,j) in Ol}(sum{c in Oc}( L[i,j] * Rc[c] *Fi[i,j,c])) + sum{(i,j) in Ol}(sum{c in Oc}(IDF * Cc[c] * L[i,j] * W[i,j,c]));
							 

#-- Restricoes --#

#Balanco de Potencia Ativa (5)
subject to BalancoPotenciaAtiva{i in Ob}:
	sum{(k,i) in Ol}(P[k,i]) - sum{(i,j) in Ol}( P[i,j] + L[i,j] * sum{c in Oc}( Rc[c]*Fi[i,j,c]) ) + PS[i]  = PD[i];

#Balanco de Potencia Reativa (6)
subject to BalancoPotenciaReativa{i in Ob}:
	sum{(k,i) in Ol}(Q[k,i]) - sum{(i,j) in Ol}( Q[i,j] + L[i,j] * sum{c in Oc}( Xc[c] *Fi[i,j,c]) ) + QS[i] = QD[i];
	
#Queda de Tensao no circuito (7)
subject to QuedaTensao{(i,j) in Ol}:
	Vsqr[i] - 2 * L[i,j] * (sum{c in Oc}(Rc[c] * Beta[i,j,c]) + sum{v in Oc}(Xc[v] * Gama[i,j,v])) - sum{n in Oc}( L[i,j]^2 * Z2[n]*  Fi[i,j,n]) - Vsqr[j] = 0;
	
#Potencia aparente (kVA)(8)
subject to PotenciaAparente{(i,j) in Ol}:
	(Vmin^2 + 0.5 * DeltaV) * Isqr[i,j] + sum{s in 1..S}(Pc[j,s]) = sum{y in 1..Y}(ms[i,j,y]*DP[i,j,y]) + sum{y in 1..Y}(ms[i,j,y]*DQ[i,j,y]);
	
# Limite de corrente (10) 
subject to LimiteCorrente{(i,j) in Ol}:
	Isqr[i,j] <= sum{c in Oc}(Icmax[c]^2 *  W[i,j,c]);

#---------------------------------------------------------------
# Limite das tensoes

subject to LimiteTensao{i in Ob}:      #(9)
	Vmin^2 <= Vsqr[i] <= Vmax^2;

subject to LinearizacaoP1_1{j in Ob}:      #(11)
	Vmin^2 + sum{s in 1..S}(xv[j,s] * DeltaV) <= Vsqr[j];

subject to LinearizacaoP1_2{j in Ob}:       #(11)
	Vsqr[j] <= Vmin^2 + DeltaV + sum{s in 1..S}(xv[j,s] * DeltaV);

#---------------------------------------------------------------
# Potencia de correcao

subject to LinearizacaoP2_1{(i,j) in Ol, s in 1..S}:         #(12)
	0 <= DeltaV * Isqr[i,j] - Pc[j,s];

subject to LinearizacaoP2_2{(i,j) in Ol, s in 1..S}:          #(12)
	DeltaV * Isqr[i,j] - Pc[j,s] <= DeltaV * sum{c in Oc}( Icmax[c]^2 *  W[i,j,c] * (1 - xv[j,s]));

subject to linearizacaoP3_1{(i,j) in Ol, s in 1..S}:          #(13)
	0 <= Pc[j,s];

subject to linearizacaoP3_2{(i,j) in Ol, s in 1..S}:           #(13)
	Pc[j,s] <= DeltaV * sum{c in Oc}( Icmax[c]^2 *  A[i,j,c,s]);

#---------------------------------------------------------------
# Variavel binaria
subject to LinearizacaoP4{j in Ob, s in 2..S}:             #(14)
	xv[j,s] <= xv[j,s-1];
	
	
#---------------------------------------------------------------	
subject to LinearizacaoP1{(i,j) in Ol}:
 Pmax[i,j] - Pmin[i,j] = P[i,j];                             #(16)
 
 
subject to LinearizacaoQ1{(i,j) in Ol}:
 Qmax[i,j] - Qmin[i,j] = Q[i,j];                            #(17)
 
 
subject to LinearizacaoP2{(i,j) in Ol}:
 Pmax[i,j] + Pmin[i,j] = sum{y in 1..Y}(DP[i,j,y]);            #(18)
 
subject to LinearizacaoQ2{(i,j) in Ol}:
 Qmax[i,j] + Qmin[i,j] =  sum{y in 1..Y}(DQ[i,j,y]);          #(19)
 
 
subject to LinearizacaoP3{(i,j) in Ol, y in 1..Y}:            #(20)
  DP[i,j,y]<= DS[i,j];
 
 subject to LinearizacaoQ3{(i,j) in Ol, y in 1..Y}:           #(21)
  DQ[i,j,y]<= DS[i,j];
 
 #------------------------------------------------------------------
 #linearização biaria
 
 #linearização fi (31)	
	 subject to LinearizacaoFi_2{(i,j) in Ol, c in Oc}:
 		0 <= -Fi[i,j,c] + Isqr[i,j]; 
 		
 	
 	#linearização fi (31)	
 subject to LinearizacaoFi{(i,j) in Ol, c in Oc}:
 	-Fi[i,j,c] + Isqr[i,j] <= M1 * (1 - W[i,j,c] );
 	
 	#linearização fi (32)
  subject to LinearizacaoFi2{(i,j) in Ol, c in Oc}:
 	Fi[i,j,c] <= M1 * W[i,j,c];
 	
  #==============================================================		
 #linearização Beta (33)
  subject to LinearizacaoBeta_2{(i,j) in Ol, c in Oc}:
 	0 <=  -Beta[i,j,c] + P[i,j];
 	
 #linearização Beta (33)
  subject to LinearizacaoBeta{(i,j) in Ol, c in Oc}:
 	-Beta[i,j,c] + P[i,j] <= M2 * (1 - W[i,j,c]);
 	
  #linearização Beta (34)
  subject to LinearizacaoBeta2{(i,j) in Ol, c in Oc}:
 	Beta[i,j,c] <= M2 *  W[i,j,c];
 	
 #==============================================================	
  
#linearização Gama (35)
  subject to LinearizacaoGama_2{(i,j) in Ol, c in Oc}:
 	0 <= -Gama[i,j,c] + Q[i,j];	
  	
 #linearização Gama (35)
  subject to LinearizacaoGama{(i,j) in Ol, c in Oc}:
 	-Gama[i,j,c] + Q[i,j] <= M3 * (1 - W[i,j,c]);	
 	
 #linearização Gama (36)
  subject to LinearizacaoGama2{(i,j) in Ol, c in Oc}:
 	Gama[i,j,c] <= M3 * W[i,j,c];	
 	
#==============================================================	
 	
 #linearização A (37)
 	subject to LinearizacaoA{(i,j) in Ol, c in Oc, s in 1..S}:
 		A[i,j,c,s] <= W[i,j,c];
 		
 #linearização A (38)
 	subject to LinearizacaoA2{(i,j) in Ol, c in Oc, s in 1..S}:
 		A[i,j,c,s] <= xv[j,s];
 	
 #linearização A (39)
 	subject to LinearizacaoA3{(i,j) in Ol, c in Oc, s in 1..S}:
 		W[i,j,c] + xv[j,s] -1 <= A[i,j,c,s];
 	
 #linearização A (39)
 	subject to LinearizacaoA3_1{(i,j) in Ol, c in Oc, s in 1..S}:
 		A[i,j,c,s] <= 1;
 	
 #==============================================================
 

 
 
 
 
 