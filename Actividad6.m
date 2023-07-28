%Seminario de Inteligencia Artificial I

%Octave 5.2.0

clear all %Limpiar datos

load data_noLineal %Cargar datos de archivo

n = size(X,1);
f1 = @(w1,w2) (1/(2*n))*sum((Y - (w1*exp(w2*X))).^2); %Funcion

hold on %Multiples datos en una gr�fica
grid on %Mostrar las lineas de la gr�fica

%Seleccion de padres
function [n] = Selection(apt)
  totalAptitud = sum(apt); %Suma de las aptitudes
  aptitudSize = numel(apt); %Tama�o de las aptitudes
  
  random = rand(); %Numero aleatorio
  selectionSum = 0; %Actual de la suma
  
  for i=1:aptitudSize
    selectionSum = selectionSum + apt(i)/totalAptitud; %Sumar la aptitud actual
    
    if selectionSum >= random
      n = i; %Aptitud seleccionada
      return;
    endif
    
  end
  
  n = aptitudSize; %Seleccionar la ultima aptitud
end

%Cruzar padres
function [y1,y2] = Cross (x1,x2)
  xSize = numel(x1); %Tama�o del padre
  crossPoint = randi([1,xSize]); %Seleccionar punto de cruce
  
  y1 = x1; %Copiar el padre
  y2 = x2; %Copiar el padre
  
  y1(crossPoint:end) = x2(crossPoint:end); %Combinar atributos desde el punto de cruce de x2
  y2(crossPoint:end) = x1(crossPoint:end); %Combinar atributos desde el punto de cruce de x1
  
end

%Mutar hijos
function [y] = Mutation (s,xl,xu)
  [row,column] = size(s); %Tama�o de los hijos
  mutationProbability = 0.01; %Probabilidad de mutaci�n
  y = s; %Copiar hijos
    
  for i=1:column
    for j=1:row
      if rand() < mutationProbability
        y(j,i) = xl(j)+(xu(j)-xl(j))*rand(); %Mutar coordenada de hijo
      endif
      
    end
    
  end
  
end

%Algoritmo Genetico

plot(X,Y,'ro','LineWidth',2,'MarkerSize',10); %Graficar valores

%Inicializacion y aptitud
xl = [0, 0]'; %Limites inferiores
xu = [1, 1]'; %Limites superiores

n = 150; %Cantidad de hijos
d = 2; %Atributos del hijo
g = 30; %Cantidad de generaciones

x = zeros(d,n); %Inicializar padres en 0
aptitud = zeros(1,n); %Inicializar aptitudes en 0

%Inicializar padres aleatoriamente
for i=1:n
  x(:,i) = xl+(xu-xl).*rand(d,1);
end

%Total de generaciones
for k=1:g
  %Calcular la aptitud de cada individuo
  for i=1:n
    functionAptitud = f1(x(1,i),x(2,i));
  
    if functionAptitud >= 0
      aptitud(i) = 1/(1+functionAptitud);
    else
      aptitud(i) = 1 + abs(functionAptitud);
    endif
  end
  
  %Conjunto vac�o de hijos
  sons = zeros(d,n);

  %Generar hijos
  for i=1:2:n
    z1 = Selection(aptitud);
    z2 = Selection(aptitud);
    
    while z1 == z2
      z2 = Selection(aptitud);
    end
    
    [son1,son2] = Cross(x(:,z1),x(:,z2)); %Genear dos hijos
    
    sons(:,i) = son1;
    sons(:,i+1) = son2;
    
  end

  sons = Mutation(sons,xl,xu); %Mutar hijos
  
  x = sons; %Los hijos son los padres
        
end

%Calcular la aptitud de cada individuo
for i=1:n
  functionAptitud = f1(x(1,i),x(2,i));

  if functionAptitud >= 0
    aptitud(i) = 1/(1+functionAptitud);
  else
    aptitud(i) = 1 + abs(functionAptitud);
  endif
end

[max,position] = max(aptitud); %Obtener la posicion

Xp = 0:0.1:10; %Valores a graficar
Yp = x(1,position)*exp(x(2,position)*Xp); %Evaluar valores

%Graficar minimo
plot(Xp,Yp,'b-','LineWidth',3,'MarkerSize',10)

%Mostrar resultados
disp("Minimo");
disp(x(1,position));
disp(x(2,position));
disp(f1(x(1,position),x(2,position)));
