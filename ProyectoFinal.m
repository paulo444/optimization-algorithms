%Seminario de Inteligencia Artificial I

%Octave 5.2.0

clear all %Limpiar datos

img = imread('Image_1.png'); %Cargar imagen a buscar
temp = imread('Template.png'); %Cargar plantilla

%Imagen y plantilla a escala de grises
img_g = rgb2gray(img);
temp_g = rgb2gray(temp);

%Correlaci�n cruzada normalizada
function val = NCC(img,temp,x,y)
    [H,W] = size(temp);
	
    sum_img = 0.0;
    sum_temp = 0.0;
    sum_2 = 0.0;
    
    x = round(x);
    y = round(y);
        
    for i=1:W
        for j=1:H
            sum_img = sum_img + double(img(y+j,x+i))^2;
            sum_temp = sum_temp + double(temp(j,i))^2;
            sum_2 = sum_2 + double(img(y+j,x+i))*double(temp(j,i));
        end
    end

    val = sum_2/(sqrt(double(sum_img))*sqrt(double(sum_temp)));
end    

%Funcion Seleccion
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

G = 15; %Numero de generaciones
N = 50; %Numero de abejas
D = 2; %Tama�o de dimension

L = 30; %Numero de intentos
Pf = 30; %Numero de abejas trabajadoras
Po = N-Pf; %Numero de abejas observadoras

[img_H,img_W] = size(img_g); %Tama�o de imagen
[temp_H,temp_W] = size(temp_g); %Tama�o de plantilla

hold on %Multiples datos en una gr�fica
grid on %Mostrar las lineas de la gr�fica

%Artificial Bee Colony

xl = [1 1]'; %Valor bajo en x
xu = [img_H-temp_H img_W-temp_W]'; %Valor alto en x

x = zeros(D,Pf); %Inicializar en cero las coordenadas
l = zeros(1,Pf); %Inicializar en cero los intentos
aptitud = zeros(1,Pf); %Inicializar en celo la aptitud
fitness = zeros(1,Pf); %Inicializar en cero las evaluaciones

%Inicializar aleatoriamente los valores
for i=1:Pf
  x(:,i) = xl+(xu-xl).*rand(D,1);
  fitness(i) = NCC(img_g,temp_g,x(2,i),x(1,i));
  
  %Calcular aptitud
  if fitness(i) >= 0
    aptitud(i) = 1/(1+fitness(i));
  else
    aptitud(i) = 1 + abs(fitness(i));
  end
end

%Generaciones
for g=1:G
  disp("GENEARACION")
  %Abejas Empleadas
  for i=1:Pf
    %Obtener fuente de alimento
    k = randi([1, Pf]);
    
    while k == i
      k = randi([1, Pf]);
    end
    
    j = randi([1, D]); %Valor aleatorio entre la dimension del problema
    phi = 2*rand()-1; %Valor aleatorio entre -1 y 1
    
    v = x(:,i); %Respaldar vector
    v(j) = x(j,i) + phi*(x(j,i)-x(j,k));
    
    %Funcion Penalizaci�n Empleadas
    for o=1:D
      if xl(o) < v(o) && v(o) < xu(o)
        %Pertence al rango
      else
        v(o) = xl(o)+(xu(o)-xl(o))*rand();
      end
    end
    
    %Evaluar en funcion
    currentFunctionV = NCC(img_g,temp_g,v(2),v(1));
    currentFunctionX = NCC(img_g,temp_g,x(2,i),x(1,i));
    
    %Asignar intentos y nuevo valor
    if currentFunctionV > currentFunctionX
      x(:,i) = v;
      fitness(i) = currentFunctionV;
      l(i) = 0;
    else
      l(i) = l(i) + 1;
    end
    
  end
  
  %Calcular aptitudes
  for i=1:Pf
    if fitness(i) >= 0
      aptitud(i) = 1/(1+fitness(i));
    else
      aptitud(i) = 1 + abs(fitness(i));
    end
  end
  
  %Abejas Observadoras
  for i=1:Po
    m = Selection(aptitud); %Abeja trabajadora
    
    %Obtener fuente de alimento
    k = randi([1, Pf]);
    
    while k == m
      k = randi([1, Pf]);
    end
    
    j = randi([1, D]); %Valor aleatorio entre la dimension del problema
    phi = 2*rand()-1; %Valor aleatorio entre -1 y 1
    
    v = x(:,m); %Respaldar vector
    v(j) = x(j,m) + phi*(x(j,m)-x(j,k));
    
    %Funcion Penalizaci�n Observadoras
    for o=1:D
      if xl(o) < v(o) && v(o) < xu(o)
        %Pertence al rango
      else
        v(o) = xl(o)+(xu(o)-xl(o))*rand();
      end
    end
    
    %Evaluar en funcion
    currentFunctionV = NCC(img_g,temp_g,v(2),v(1));
    currentFunctionX = NCC(img_g,temp_g,x(2,m),x(1,m));
    
    %Asignar intentos y nuevo valor
    if currentFunctionV > currentFunctionX
      x(:,m) = v;
      fitness(m) = currentFunctionV;
      l(m) = 0;
    else
      l(m) = l(m) + 1;
    end
    
  end
  
  %Abejas Exploradoras
  for i=1:Pf
    %Buscar otra posicion si alcanzo el maximo de intentos
    if l(i) > L
      fitness(i) = NCC(img_g,temp_g,x(2,i),x(1,i));
      x(:,i) = xl+(xu-xl).*rand(D,1);
      l(i) = 0;
    end
    
  end
  
end

[~,bestPosition] = max(fitness); %Obtener mejor posicion
xPosition = round(x(2,bestPosition));
yPosition = round(x(1,bestPosition));

%Graficar minimo
imshow(img)

line([xPosition xPosition+temp_W], [yPosition yPosition],'Color','g','LineWidth',3);
line([xPosition xPosition], [yPosition yPosition+temp_H],'Color','g','LineWidth',3);
line([xPosition+temp_W xPosition+temp_W], [yPosition yPosition+temp_H],'Color','g','LineWidth',3);
line([xPosition xPosition+temp_W], [yPosition+temp_H yPosition+temp_H],'Color','g','LineWidth',3);

%Mostrar resultados
disp("Minimo");
disp(xPosition);
disp(yPosition);
disp(NCC(img_g,temp_g,xPosition,yPosition));