# SIFIN Bolivia — Sistema Inteligente de Planificación y Distribución Financiera

Proyecto académico de **Métodos Numéricos** para la mención de **Ingeniería de Sistemas**. Implementa una aplicación web interactiva que modela la distribución de recursos entre cuatro unidades de una empresa boliviana y compara los métodos de Gradientes Conjugados (GC) y Gradientes Conjugados Precondicionado (GCP).

> **Aviso académico:** las matrices, coeficientes y necesidades financieras incluidos son datos simulados y construidos con fines didácticos. No representan presupuestos oficiales ni información de una empresa real. Los resultados no constituyen asesoramiento financiero.

## 1. Requisitos del desafío cubiertos

- Sistema lineal cuadrado de dimensión 4: `A x = b`.
- Matriz simétrica y definida positiva, validada mediante factorización de Cholesky.
- Tres escenarios coherentes con las mismas cuatro variables: Ideal, Bajo estrés financiero y Mal condicionado.
- Solución numérica de referencia mediante Cholesky.
- Implementación de GC simple y GCP con precondicionador diagonal de Jacobi.
- Vector inicial editable, tolerancia configurable y máximo de iteraciones.
- Criterio de parada: `||r_k||₂ / ||b||₂ < tolerancia` (por defecto `10⁻⁶`).
- Errores absoluto, relativo y porcentual por variable; residuo relativo global.
- Número de condición de `A` y de la matriz simétrica equivalente precondicionada.
- Comparación de iteraciones, residuos y convergencia.
- Entrada dinámica de matriz `A` y vector `b`.
- Módulo QA con cuatro pruebas requeridas.
- Conclusiones, interpretación financiera y referencias.

## 2. Tecnologías

- HTML5: estructura y contenido.
- CSS3: diseño adaptable (responsive).
- JavaScript: interfaz, validación y algoritmos numéricos.
- No requiere Node.js, base de datos, instalación de paquetes ni backend para ejecutar la demostración.

La tipografía se carga desde Google Fonts cuando hay conexión a Internet. Si no hay conexión, el navegador utiliza tipografías de respaldo.

## 3. Estructura del proyecto

```text
SIFIN-Bolivia/
├── index.html
├── README.md
├── css/
│   └── styles.css
└── js/
    ├── math.js
    └── app.js
```

### Responsabilidad de los archivos

- `index.html`: estructura de los diez módulos de la página.
- `css/styles.css`: identidad visual, componentes y adaptación a móviles.
- `js/math.js`: operaciones vectoriales y matriciales, Cholesky, autovalores por Jacobi, número de condición, validación, GC y GCP.
- `js/app.js`: escenarios, edición de datos, navegación, presentación de resultados, errores y pruebas QA.

## 4. Instalación y ejecución local

### Opción A: Visual Studio Code

1. Descarga o clona el repositorio.
2. Descomprime la carpeta si descargaste un ZIP.
3. Abre Visual Studio Code.
4. Selecciona **Archivo → Abrir carpeta** y elige `SIFIN-Bolivia`.
5. Comprueba que existan `index.html`, `css/styles.css`, `js/math.js` y `js/app.js`.
6. Abre `index.html` en el navegador. También puedes utilizar la extensión **Live Server** de Visual Studio Code y seleccionar **Open with Live Server**.

### Opción B: navegador

Abre directamente el archivo `index.html` con Chrome, Firefox o Edge. La aplicación funciona en el navegador sin servidor.

## 5. Arquitectura y flujo de ejecución

La aplicación está organizada en tres capas sencillas:

1. **Presentación:** `index.html` y `css/styles.css` muestran los módulos, formularios, tablas y diagnósticos.
2. **Aplicación:** `js/app.js` gestiona los escenarios, datos editables, navegación, ejecución y renderizado de resultados.
3. **Núcleo numérico:** `js/math.js` realiza los cálculos sin depender de librerías externas.

Flujo general:

```text
Seleccionar escenario / editar A, b y x₀
                  ↓
       Validar dimensiones y datos
                  ↓
       Comprobar simetría y SPD
                  ↓
   Solución de referencia (Cholesky)
                  ↓
          ┌───────┴───────┐
          ↓               ↓
      GC simple       GCP + Jacobi
          ↓               ↓
          └───────┬───────┘
                  ↓
    Errores, residuos y condicionamiento
                  ↓
       Tablas, diagnósticos y QA
```

## 6. Modelo financiero

Se considera una empresa boliviana simulada con cuatro unidades:

- `x₁`: Producción.
- `x₂`: Logística y transporte.
- `x₃`: Comercialización.
- `x₄`: Administración.

El sistema `A x = b` es una representación lineal simplificada de la coordinación y asignación de recursos. Los términos diagonales representan el peso o capacidad propia de cada unidad; los términos fuera de la diagonal representan interdependencias entre unidades. El vector `b` representa necesidades financieras estimadas para el periodo.

Los coeficientes son supuestos didácticos. El modelo no debe interpretarse como una contabilidad real ni como un modelo de optimización financiera con restricciones empresariales completas.

### Escenarios

- **Ideal:** operación estable y matriz moderadamente condicionada.
- **Bajo estrés financiero:** se elevan las necesidades simuladas de producción y logística y se modifican algunos acoplamientos.
- **Mal condicionado:** se utilizan escalas diagonales muy diferentes para ilustrar el mal escalamiento. Se aplica Jacobi como precondicionador diagonal.

Los escenarios mantienen las mismas cuatro variables y la misma estructura 4 × 4.

## 7. Métodos numéricos

### Solución de referencia

Se utiliza factorización de Cholesky `A = L Lᵀ`, válida para matrices simétricas definidas positivas. Luego se resuelven los sistemas triangulares `L y = b` y `Lᵀ x = y`.

### Gradientes Conjugados (GC)

Se implementa el algoritmo iterativo estándar con residuo inicial `r₀ = b − A x₀`, dirección inicial `p₀ = r₀` y actualización de solución, residuo y dirección conjugada.

### Gradientes Conjugados Precondicionado (GCP)

Se utiliza el precondicionador diagonal de Jacobi `M = diag(A)`. En cada iteración se resuelve `M zₖ = rₖ`, equivalente a dividir cada componente del residuo por el correspondiente elemento diagonal. Se actualizan `αₖ`, `βₖ`, `xₖ`, `rₖ` y `pₖ`.

### Número de condición precondicionado

Para informar el condicionamiento se calcula el número de condición espectral de la matriz simétrica equivalente `M⁻¹/² A M⁻¹/²`, que tiene los mismos autovalores que `M⁻¹ A` cuando `M` es definida positiva. Los autovalores se aproximan mediante rotaciones de Jacobi.

### Errores

Para cada componente:

- Error absoluto: `Ea,i = |xref,i − xGCP,i|`.
- Error relativo: `Er,i = Ea,i / |xref,i|`.
- Error porcentual: `Ep,i = Er,i × 100`.
- Si `xref,i = 0`, el error relativo y porcentual de esa componente se informa como no definido.

El residuo relativo global es `||b − A xGCP||₂ / ||b||₂`.

## 8. Módulos de la página

1. Inicio.
2. Caso de estudio.
3. Formulación matemática.
4. Entrada de datos.
5. Tres escenarios.
6. Solución de referencia y comparación GC vs GCP.
7. Análisis de errores.
8. Módulo de pruebas QA.
9. Conclusiones y referencias.

La navegación contiene diez módulos independientes: Inicio, Caso de estudio, Formulación, Entrada de datos, Escenarios, Solución de referencia, GC vs GCP, Análisis de errores, QA y Conclusiones/Referencias.

## 9. Pruebas QA incluidas

El módulo QA ejecuta estas verificaciones:

1. **Carga de escenarios:** valida simetría, definición positiva y resolución directa para Ideal, Estrés y Mal condicionado.
2. **Matriz no válida:** altera intencionalmente un coeficiente para producir una matriz no simétrica y comprueba que la validación capture el error.
3. **Edición dinámica:** modifica un componente del vector `b` y verifica que la solución de referencia cambie al recalcular.
4. **Comparación de precondicionamiento:** comprueba en el escenario mal condicionado que GCP converja y use un número de iteraciones menor o igual que GC. Si el resultado no se cumple, la prueba se marca como no aprobada.

Las pruebas se ejecutan desde el botón **Ejecutar todas las pruebas** de la vista QA.

## 10. Publicación y repositorio

### GitHub

1. Inicia sesión en GitHub y crea un repositorio público llamado `SIFIN-Bolivia`.
2. Sube el contenido de esta carpeta conservando la estructura de directorios.
3. Incluye commits descriptivos durante el desarrollo.
4. Para publicar la página, entra a **Settings → Pages**.
5. En **Build and deployment**, selecciona **Deploy from a branch**.
6. Selecciona la rama `main` y la carpeta `/ (root)`; guarda los cambios.
7. Espera a que GitHub Pages publique el sitio. La URL normalmente tendrá la forma `https://TU-USUARIO.github.io/SIFIN-Bolivia/`.
8. Comprueba el enlace en una ventana privada o en otro dispositivo.

No se incluye una URL pública ficticia: el enlace final se genera cuando el proyecto se publica desde tu propia cuenta.

## 11. Referencias

- Suñagua S., P. (2020). *Método de Gradientes Conjugados Precondicionado*. Revista Boliviana de Matemática – UMSA, 04, 2–7. Registro institucional: https://iimat.umsa.bo/publications/
- Golub, G. H., & Van Loan, C. F. (2013). *Matrix Computations* (4th ed.). Johns Hopkins University Press.
- Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems* (2nd ed.). SIAM.

## 12. Autoría y alcance

Proyecto académico individual — Mención: Ingeniería de Sistemas — Universidad Mayor de San Andrés. Completa en la portada los datos personales que tu docente solicite antes de la presentación. Los cálculos se ejecutan localmente en el navegador; no se recopilan ni transmiten datos.
