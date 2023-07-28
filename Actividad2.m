%Seminario de Inteligencia Artificial I

%Octave 5.2.0

clear all %Limpiar datos

functionNumber = input("Funcion a optimizar: "); %Numero de funcion a ejecutar
methodNumber = input("Metodo a utilizar: "); %Metodo a utilizar

f1 = @(x,y) x.*e.^(-x.^2-y.^2); %Funcion 1

f2 = @(x,y) (x-2).^2 + (y-2).^2; %Funcion 2

hold on %Multiples datos en una gr�fica
grid on %Mostrar las lineas de la gr�fica

% Funcion de recombinaci�n
function y = Recombination (x1,x2)
    y = 0.5*(x1+x2);
end

G = 100; %Numero de generaciones
mu = 50; %Numero de padres
D = 2; %Tama�o de valores
lambda = 60; %Numero de hijos

%Seleccionar funci�n
switch functionNumber
  %Funci�n 1
  case 1
    xl = -2; %Parte baja en x
    xu = 2; %Parte alta en x
    yl = -2; %Parte baja en y
    yu =2; %Parte alta en y
    
    [x,y] = meshgrid(xl:0.2:xu,yl:0.2:yu); %Coordenadas a graficar
    z = f1(x,y); %Valores de la funcion evaluada

    contour(x,y,z) %Graficar
    
    xl = [-2 -2]'; %Parte baja
    xu = [2 2]'; %Parte alta
    
    x = zeros(D,mu+lambda); %Inicializar en ceros los padres
    sigma = zeros(D,mu+lambda); %Inicializar en ceros los sigmas
    fitness = zeros(1,mu+lambda); %Inicializar en ceros el fitness
    
    %Inicializar aleatoriamente
    for i=1:mu
      x(:,i) = xl+(xu-xl).*rand(D,1);
      sigma(:,i) = 0.5*rand(D,1);
      fitness(i) = f1(x(1,i),x(2,i));
    end
    
    %Generaciones
    for k=1:G
      %Crear hijos
      for i=1:lambda
        %Elegir padres
        r1 = randi([1 mu]);
        r2 = randi([1 mu]);
        while r1==r2
            r2 = randi([1 mu]);
        end
        
        %Recombinar padres
        x(:,mu+i) = Recombination(x(:,r1),x(:,r2));
        sigma(:,mu+i) = Recombination(sigma(:,r1),sigma(:,r2));
        
        %Calcular normalizaci�n
        r = normrnd(0,sigma(:,mu+i),[D 1]);
        x(:,mu+i) = x(:,mu+i) + r;
      
        %Evaluar fitness
        fitness(mu+i) = f1(x(1,mu+i),x(2,mu+i));
        
      end
      
      switch methodNumber
        %(? + ?)-ES
        case 1
          %Ordenar hijos y padres
          [~,I] = sort(fitness);
          x = x(:,I);
          sigma = sigma(:,I);
          fitness = fitness(I);
        
        %(?,?)-ES
        case 2
          %Inicializar hijos
          y = zeros(D,lambda);
          ySigma = zeros(D,lambda);
          yFitness = zeros(1,lambda);
          
          %Obtener hijos
          for j=1:lambda
            y(:,j) = x(:,mu+j);
            ySigma(:,j) = sigma(:,mu+j);
            yFitness(:,j) = fitness(1,mu+j);
          end
          
          %Ordenar hijos
          [~,I] = sort(yFitness);
          y = y(:,I);
          ySigma = ySigma(:,I);
          yFitness = yFitness(:,I);
          
          %Copiar hijos ordenados
          for j=1:lambda
            x(:,j) = y(:,j);
            sigma(:,j) = ySigma(:,j);
            fitness(:,j) = yFitness(:,j);
          end
          
          %Borar padres
          for j=1:lambda
            x(:,mu+j) = zeros(D,1);
            sigma(:,mu+j) = zeros(D,1);
            fitness(:,mu+j) = zeros(1,1);
          end

      endswitch
      
    end
    
    %Graficar minimo
    plot(x(1,1),x(2,1),'rx','LineWidth',2,'markerSize',10) %Mostrar el optimo
    
    %Mostrar resultados
    disp("Minimo");
    disp(x(1,1));
    disp(x(2,1));
    disp(f1(x(1,1),x(2,1)));
  
  %Funci�n 2
  case 2
    xl = -5; %Parte baja en x
    xu = 5; %Parte alta en x
    yl = -5; %Parte baja en y
    yu =5; %Parte alta en y
    
    [x,y] = meshgrid(xl:0.2:xu,yl:0.2:yu); %Coordenadas a graficar
    z = f2(x,y); %Valores de la funcion evaluada

    contour(x,y,z) %Graficar
    
    xl = [-5 -5]'; %Parte baja
    xu = [5 5]'; %Parte alta
    
    x = zeros(D,mu+lambda); %Inicializar en ceros los padres
    sigma = zeros(D,mu+lambda); %Inicializar en ceros los sigmas
    fitness = zeros(1,mu+lambda); %Inicializar en ceros el fitness
    
    %Inicializar aleatoriamente
    for i=1:mu
      x(:,i) = xl+(xu-xl).*rand(D,1);
      sigma(:,i) = 0.5*rand(D,1);
      fitness(i) = f2(x(1,i),x(2,i));
    end
    
    %Generaciones
    for k=1:G
      %Crear hijos
      for i=1:lambda
        %Elegir padres
        r1 = randi([1 mu]);
        r2 = randi([1 mu]);
        while r1==r2
            r2 = randi([1 mu]);
        end
        
        %Recombinar padres
        x(:,mu+i) = Recombination(x(:,r1),x(:,r2));
        sigma(:,mu+i) = Recombination(sigma(:,r1),sigma(:,r2));
        
        %Calcular normalizaci�n
        r = normrnd(0,sigma(:,mu+i),[D 1]);
        x(:,mu+i) = x(:,mu+i) + r;
      
        %Evaluar fitness
        fitness(mu+i) = f2(x(1,mu+i),x(2,mu+i));
        
      end
      
      switch methodNumber
        %(? + ?)-ES
        case 1
          %Ordenar hijos y padres
          [~,I] = sort(fitness);
          x = x(:,I);
          sigma = sigma(:,I);
          fitness = fitness(I);
        
        %(?,?)-ES
        case 2
          %Inicializar hijos
          y = zeros(D,lambda);
          ySigma = zeros(D,lambda);
          yFitness = zeros(1,lambda);
          
          %Obtener hijos
          for j=1:lambda
            y(:,j) = x(:,mu+j);
            ySigma(:,j) = sigma(:,mu+j);
            yFitness(:,j) = fitness(1,mu+j);
          end
          
          %Ordenar hijos
          [~,I] = sort(yFitness);
          y = y(:,I);
          ySigma = ySigma(:,I);
          yFitness = yFitness(:,I);
          
          %Copiar hijos ordenados
          for j=1:lambda
            x(:,j) = y(:,j);
            sigma(:,j) = ySigma(:,j);
            fitness(:,j) = yFitness(:,j);
          end
          
          %Borar padres
          for j=1:lambda
            x(:,mu+j) = zeros(D,1);
            sigma(:,mu+j) = zeros(D,1);
            fitness(:,mu+j) = zeros(1,1);
          end

      endswitch
      
    end
    
    %Graficar minimo
    plot(x(1,1),x(2,1),'rx','LineWidth',2,'markerSize',10) %Mostrar el optimo
    
    %Mostrar resultados
    disp("Minimo");
    disp(x(1,1));
    disp(x(2,1));
    disp(f2(x(1,1),x(2,1)));

  otherwise
    disp("Numero no valido")

endswitch