# Modelo de Colapso Esférico (SCM)

Este proyecto implementa un programa para calcular la evolución de las perturbaciones esféricas de densidad mediante el modelo de colapso esférico (SCM). Ofrece herramientas para explorar diferentes modelos cosmológicos y estudiar su impacto en la formación de estructuras a gran escala.

## Características principales

- **Modelos cosmológicos:** Los modelos se definen en el archivo `lib/DE-models.cpp`.
  - **Modelos implementados:**
    - Modelos EdS.
    - Modelos LCDM.
    - Modelo CPL.
    - Modelo CPLF (CPL Phanton).
    - Modelo neDE: Energía oscura temprana.
    - Modelo 2EXP: doble exponecial.
    - Modelo con ecuación de estado dinámica (numérica) W. "Modelo Melissa"; devgf.
    - Modelo con ecuación de estado dinámica (numérica) W. "Non Abelian"; nonabelian.
    - Modelo con ecuación de estado dinámica (numérica) W. "Non Abelian"; twoform.
- **Ecuaciones del modelo SCM:** Implementadas en `lib/scm.cpp`.
- **Selección de modelos:** Configurable mediante `modelSelect.cpp`.
- **Rutinas principales:** El archivo `main.cpp` contiene las rutinas para calcular las cantidades relacionadas con el modelo SCM.

## Estructura del proyecto

```
SCM-v1.55/
├── lib/
│   ├── DE-models.cpp     # Definición de modelos cosmológicos
│   ├── scm.cpp           # Implementación de ecuaciones del SCM
│   └── modelSelect.cpp   # Selección del modelo a analizar
├── main.cpp              # Rutinas principales para el cálculo
├── include/              # Archivos de cabecera
├── output/               # Resultados generados
└── readme.md             # Documentación original
```

## Instalación

1. Clona el repositorio o extrae el contenido del archivo ZIP.
2. Asegúrate de tener un compilador C++ compatible (por ejemplo, `g++`).
3. Compila el proyecto:

    ```bash
    make comp
    ```
    o

   ```bash
   g++ -o scm main.cpp lib/*.cpp -Iinclude
   ```

## Uso

1. Edita `modelSelect.cpp` para seleccionar el modelo cosmológico deseado.
2. Ejecuta el programa:

    ```bash
    make run
    ```

3. Los resultados se guardarán en la carpeta `output/` para su análisis.

## Limpiar archivos

Para limpiar los archivos de compilación usar:

    
    make clean
    

## Personalización

- **Nuevos modelos cosmológicos:** Puedes agregar nuevos modelos en `lib/DE-models.cpp` siguiendo el formato existente.
- **Ecuaciones adicionales:** Para modificar o extender las ecuaciones del modelo SCM, edita `lib/scm.cpp`.

## Contribuciones

Las contribuciones al proyecto son bienvenidas. Si tienes sugerencias o mejoras, abre un issue o envía un pull request.

## Licencia

Este proyecto está bajo una licencia abierta. Consulta el archivo `LICENSE` para más detalles.


