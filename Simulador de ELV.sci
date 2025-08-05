clc
clear
clf

// Simulador de Equilíbrio Líquido-Vapor Binário com UNIFAC, Wilson e Virial
// 
// Este programa realiza a simulação do equilíbrio líquido-vapor (ELV) para
// uma mistura binária utilizando modelos termodinâmicos avançados. São aplicados
// dois modelos de coeficiente de atividade — UNIFAC (baseado em grupos funcionais)
// e Wilson (baseado em parâmetros de interação binária) — para estimar a não
// idealidade da fase líquida. Além disso, o modelo de virial truncado é empregado
// para corrigir a não idealidade da fase vapor via coeficientes de fugacidade.
// O algoritmo ajusta iterativamente a temperatura de equilíbrio para diferentes
// composições líquidas, comparando os resultados com dados experimentais e
// fornecendo gráficos e análise de erro.


//Parâmtros de Antoine
A = [12.2917,3803.93,-41.68;
    9.1325,2766.63,-50.50]
    
//Propriedades termodinâmicas (w,Tc,Pc)
W = [0.635,516.2,63.83;
    0.213,553.4,40.73]
    
T = 310
P = 1.013
R = 8.314

//Matrizes do modelo UNIFAC - (Tiradas da tabela em ANEXOS)

Qi = [2.588,3.24]; Ri = [2.575,4.044]

G = [0,0,986.5;
    0,0,986.5;
    156.4,156.4,0]
    
E = [0.327666151,0;
    0.208655332,1;
    0.463678516,0]

//Função para cálculo de gamma

//Tebela dos artigos

U = [348.85 0.9973 0.8973
    347.55 0.9890 0.8431
    346.55 0.9742 0.7960
    344.75 0.9499 0.7230
    342.75 0.9119 0.6520
    341.45 0.8775 0.5942
    340.05 0.8280 0.5361
    339.35 0.7849 0.5084
    338.55 0.6917 0.4763
    338.35 0.6425 0.4742
    337.95 0.5693 0.4668
    337.85 0.5401 0.4564
    337.65 0.4345 0.4482
    337.75 0.3227 0.4457
    337.85 0.2734 0.4530
    338.75 0.1393 0.4680
    340.55 0.0648 0.4281
    342.15 0.0526 0.4052
    345.75 0.0352 0.3595
    347.65 0.0303 0.3160
    348.65 0.0245 0.2870
    350.15 0.0192 0.2332
    351.35 0.0157 0.2066
    352.75 0.0144 0.1341
    353.45 0.0102 0.0688]

function f=gam(T, G, E, Qi, Ri, x)
    for n=1:3
        for k=1:3
            t(k,n) = exp(-G(n,k)/T)
        end,
    end,
    
    B = t*E
    
    dl(1) = x*Qi(1)*E(1,1)/(x*Qi(1)+(1-x)*Qi(2))+(1-x)*Qi(2)*E(1,2)/(x*Qi(1)+(1-x)*Qi(2))
    dl(2) = x*Qi(1)*E(2,1)/(x*Qi(1)+(1-x)*Qi(2))+(1-x)*Qi(2)*E(2,2)/(x*Qi(1)+(1-x)*Qi(2))
    dl(3) = x*Qi(1)*E(3,1)/(x*Qi(1)+(1-x)*Qi(2))+(1-x)*Qi(2)*E(3,2)/(x*Qi(1)+(1-x)*Qi(2))
    
    S = t*dl
    
    J(1) = Ri(1)/(Ri(1)*x+Ri(2)*(1-x)); L(1) = Qi(1)/(Qi(1)*x+Qi(2)*(1-x))
    J(2) = Ri(2)/(Ri(1)*x+Ri(2)*(1-x)); L(2) = Qi(2)/(Qi(1)*x+Qi(2)*(1-x))
    
    lnyc(1) = 1-J(1)+log(J(1))-5*Qi(1)*(1-J(1)/L(1)+log(J(1)/L(1)))
    lnyc(2) = 1-J(2)+log(J(2))-5*Qi(2)*(1-J(2)/L(2)+log(J(2)/L(2)))
    
    for i=1:3
        k(1,i) = dl(i)*B(i,1)/S(i)-E(i,1)*log(B(i,1)/S(i))
        k(2,i) = dl(i)*B(i,2)/S(i)-E(i,2)*log(B(i,2)/S(i))
    end
    
    lnyr(1) = Qi(1)*(1-(k(1,1)+k(1,2)+k(1,3)))
    lnyr(2) = Qi(2)*(1-(k(2,1)+k(2,2)+k(2,3)))
    
    f(1,1) = exp(lnyr(1)+lnyc(1)),f(1,2) = exp(lnyr(2)+lnyc(2))
endfunction

function f=wilson(T, V1, V2, x1, x2, a12, a21, R)
    A12 = V2/V1*exp(-a12/(R*T)), A21 = V1/V2*exp(-a21/(R*T))
    A11 = 1, A22 = 1
    x2 = 1-x1
    lny1 = 1-log(A11*x1+A12*x2)-(A11*x1/(A11*x1+A12*x2)+A21*x2/(A21*x1+A22*x2))
    lny2 = 1-log(A21*x1+A22*x2)-(A12*x1/(A11*x1+A12*x2)+A22*x2/(A21*x1+A22*x2))
    f(1) = exp(lny1), f(2) = exp(lny2)
    
endfunction

//Função para virial truncada no 2 termo

function f=virial(y1, y2, T, P, W, R)
    for i=1:2
        w = W(i,1)
        Tc = W(i,2)
        Pc = W(i,3)
        Tr = T/Tc
        B0 = 0.083-0.422/(Tr^1.6); B1 = 0.139-0.172/(Tr^4.2);
        B = (R*Tc/Pc)*(B0+w*B1)
        ANS(1,i) = B
    end
    B12 = (ANS(1,1)+ANS(1,2))/2
    Bk = 2*y1*y2*B12
    lnO1 = (P/(R*T))*(2*(y2*B12)-Bk); lnO2 = (P/(R*T))*(2*(y1*B12)-Bk);
    f(1) = exp(lnO1); f(2) = exp(lnO2);
endfunction

//Função para cálculo de Psat e Phisat por Antoine

function f=antoine(T, A, W, R)
    for i=1:2
        
        a = A(i,1); b = A(i,2); c = A(i,3);
        
        Psat = exp(a-(b/(T+c))) //Psat  = 10^(a-(b/(T+c)))
        f(i,1) = Psat
        
        w = W(i,1)
        Tc = W(i,2)
        Pc = W(i,3)
        Tr = T/Tc
        B0 = 0.083-0.422/(Tr^1.6); B1 = 0.139-0.172/(Tr^4.2);
        B = (R*Tc/Pc)*(B0+w*B1)
        lnOsat = B*Psat/(R*T)
        f(i,2) = exp(lnOsat)
        
    end
endfunction

function f=algoritmo(T, P, A, W, G, E, Qi, Ri, R, nk, U, C1, C2)
    
    T_ = T //Temperatura auxiliar
    nj = 1 //Iterações
    
    for x = 0:(1/nk):1
        
        //Condições e valores auxiliares
        y_correto = %F, n2 = 1 
        maior = 5, menor = 5, c = 1.1, T = T_
        
        if C1 == %T then //Para valores pré definidos de x
            for k = 1:nk
                x = U(k,2)
            end
        end,
        
        while y_correto == %F
            
            mudou = %T, n = 1
            
            ANT = antoine(T,A,W,R), 
            
            if C2 == %F then //Para UNIFAC
                GAM = gam(T,G,E,Qi,Ri,x)
            else //Para Wilson
                V1 = 59.08, V2 = 108.75
                a12 = 8072.3, a21 = 2098.8 
                GAM = wilson(T,V1,V2,x,(1-x),a12,a21,R)
            end,
            
            P1sat = ANT(1,1), P2sat = ANT(2,1)
            O1sat = ANT(1,2), O2sat = ANT(2,2)
            gam1 = GAM(1), gam2 = GAM(2), O1 = 1, O2 = 1
            y1 = x*P1sat*O1sat*gam1/(P*O1)
            y2 = (1-x)*P2sat*O2sat*gam2/(P*O2)
            
            while mudou == %T
                if n == 1 then
                    y = y1+y2
                    y_ = y
                end,
                
                VIR = virial(y1,y2,T,P,W,R), O1 = VIR(1), O2 = VIR(2)
                y1 = x*P1sat*O1sat*gam1/(P*O1)
                y2 = (1-x)*P2sat*O2sat*gam2/(P*O2)
                y = y1+y2
                
                if abs(y-y_) < 1E-10 then
                    mudou = %F
                end,
                
                n = n+1
                y_ = y
                
                if n> 1000 then
                    mudou = %F
                    disp("erro1")
                end,
            end,
            
            if abs(y-1) > 1E-4 then
                
                if y < 1 then
                    T = T + maior
                    menor = menor/c
                else
                    T = T - menor
                    maior = maior/c
                end,
                
                n2 = n2+1
            else
                y_correto = %T
            end,
            
            if n2> 10000 then
                y_correto = %T
                disp("erro2")
            end,
        end,
        f(1,nj) = x, f(2,nj) = T, f(3,nj) = y1, f(4,nj) = y2
        f(5,nj) = gam1, f(6,nj) = gam2
        nj = nj+1
    end
endfunction

ALG = algoritmo(T,P,A,W,G,E,Qi,Ri,R,100,U,%F,%F)
ALGW = algoritmo(T,P,A,W,G,E,Qi,Ri,R,100,U,%F,%T)

//Plot dos gráficos

X1 = ALG(1,:), Y1 = ALG(2,:), X1W = ALGW(1,:), Y1W = ALGW(2,:)
show_window(0)
plot(X1,Y1,"-k")
plot(X1W,Y1W,"-r")
legend("UNIFAC","Wilson")
xtitle("T vs. x1")
a = gca()
a.font_size = 2
a.x_label.text = "$x1$", a.x_label.font_size = 4
a.y_label.text = "$T (K)$", a.y_label.font_size = 4

Y2 = ALG(3,:), Y3 = ALG(4,:), Y2W = ALGW(3,:), Y3W = ALGW(4,:)
show_window(1)
plot(X1,Y2,"-k")
plot(X1,Y3,"-b")
plot(X1W,Y2W,"-r")
plot(X1W,Y3W,"-g")
legend("UNIFAC y1","UNIFAC y2", "Wilson y1", "Wilson y2")
xtitle("y1,y2 vs. x")
a = gca()
a.font_size = 2
a.x_label.text = "$x$", a.x_label.font_size = 4
a.y_label.text = "$y1,y2$", a.y_label.font_size = 4

show_window(2)
plot(Y2,Y1,"-k")
plot(X1,Y1,"-b")
plot(Y2W,Y1W,"-r")
plot(X1W,Y1W,"-g")
legend("UNIFAC y","UNIFAC x","Wilson y","Wilson x")
xtitle("T vs. x1 e y1")
a = gca()
a.font_size = 2
a.x_label.text = "$x1, y1$", a.x_label.font_size = 4
a.y_label.text = "$T(K)$", a.y_label.font_size = 4
a.data_bounds = [0,335;1,360]

Y4 = ALG(5,:), Y5 = ALG(6,:), Y4W = ALGW(5,:), Y5W = ALGW(6,:)
show_window(3)
plot(X1,log(Y4),"-k")
plot(X1,log(Y5),"-b")
plot(X1W,log(Y4W),"-r")
plot(X1W,log(Y5W),"-g")
legend("UNIFAC γ1", "UNIFAC γ2", "Wilson γ1", "Wilson γ2")
xtitle("ln γ vs. x1")
a = gca()
a.font_size = 2
a.x_label.text = "$x1$", a.x_label.font_size = 4
a.y_label.text = "$lnγ$", a.y_label.font_size = 4

//Função para calcular o erro

ALG_E = algoritmo(T,P,A,W,G,E,Qi,Ri,R,25,U,%T,%F)
ALG_EW = algoritmo(T,P,A,W,G,E,Qi,Ri,R,25,U,%T,%T)
X6 = U(:,2), Y6 = ALG_E(2,:), Y_E = U(:,1),Y6W = ALG_EW(2,:)
for i=1:25
    ERRO(i) = abs((Y_E(i)-Y6(i))/Y_E(i)*100)
    ERROW(i) = abs((Y_E(i)-Y6W(i))/Y_E(i)*100)
end
show_window(4)
plot(X6,ERRO,"ok")
plot(X6,ERROW,"or")
xtitle("Erro vs x1")
legend("UNIFAC","Wilson")
a = gca()
a.font_size = 2
a.x_label.text = "$x1$", a.x_label.font_size = 4
a.y_label.text = "$Erro$", a.y_label.font_size = 4


