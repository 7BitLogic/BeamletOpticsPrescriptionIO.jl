# BeamletOpticsZMX.jl

**BeamletOpticsZMX.jl** ist ein eigenständiges Import- und Konvertierungstool, mit dem optische Systeme aus **Zemax OpticStudio** (`.zmx`-Dateien) eingelesen und als vollwertige 3D-Komponentensysteme für [BeamletOptics.jl](https://github.com/7BitLogic/BeamletOptics.jl) (BMO) rekonstruiert werden können.

Das Tool basiert konzeptionell auf den Parsing- und Sequentiell-zu-Element-Gruppierungsmechanismen des Python-Projekts [`ray-optics`](https://github.com/mjhoptics/ray-optics) von Michael J. Hayford, transformiert die 2D-Flächenmodelle jedoch in die dreidimensionalen, volumetrischen Signed Distance Function (SDF) Festkörperstrukturen von `BeamletOptics.jl`.

---

## Kernfunktionen

- **Robuster ZMX-Parser**: Liest Zemax-Dateien in UTF-16LE (mit/ohne BOM), UTF-8 und ISO-8859-1. Unterstützt Einheitenumrechnung (`MM`, `IN`, `CM`, `M`), Wellenlängen, Krümmungen, Dicken, Asphärenkoeffizienten und Blenden.
- **Sequentiell-zu-Element-Assemblierung**: Gruppiert kontinuierliche Flächenketten automatisch in konkrete Baugruppen:
  - **Einzellinsen** (`SphericalLens` bzw. asphärische `Lens` mit `EvenAsphericalSurface`)
  - **Verkittete Dubletts** (`SphericalDoubletLens`)
  - **Verkittete Tripletts** (`SphericalTripletLens`)
  - **Spiegel** (`RoundPlanoMirror`, `SphericalMirror`)
  - **Aperturblenden** (`Stop`)
  - **Bildsensoren / Detektoren** (`Detector`)
- **Koordinatentransformation & Kinematik**:
  - Konvertiert die Zemax-Ausbreitungsachse (+Z) auf die optische Hauptachse von BeamletOptics (+Y).
  - Skaliert alle Dimensionen präzise in SI-Einheiten (Meter).
  - Positioniert alle Komponenten an ihren exakten kumulativen Scheitelpunktkoordinaten.
- **Glas- & Dispersionsdatenbank**:
  - Nahtlose Anbindung an `RefractiveIndex.jl` für über 1000 optische Kataloggläser (Schott, Ohara, Hoya, CDGM, Hikari, Sumita).
  - Berechnung der physikalischen Dispersion für Modellgläser ($n_d, V_d$) via Cauchy-Approximation über die Fraunhofer-Linien (F, d, C).
- **Zwei Nutzungsmodi**:
  1. **Direkter Julia-Import**: Liefert direkt ein simulierbares `BMO.System` und `Detector`.
  2. **Code-Generator & CLI**: Erzeugt lesbaren, eigenständigen Julia-Quelltext (`.jl`), der ohne Abhängigkeit zum ZMX-Importer ausgeführt werden kann.

---

## Installation & Einrichtung

```julia
using Pkg
Pkg.activate("path/to/BMO_ZMX_import")
Pkg.instantiate()
```

---

## Verwendung

### 1. Direkter Import in Julia

```julia
using BeamletOptics
using BeamletOpticsZMX

# ZMX-Datei importieren
res = import_zmx("test_data/dan_reiley/photographic-lenses-prime/US00583336-2-scaled.zmx")

# Zugriff auf das generierte BMO-System und Detektor
system = res.system
detector = res.detector

# Strahlverfolgung durchführen
source = CollimatedSource([0.0, -10.0e-3, 0.0], [0.0, 1.0, 0.0], 5.0e-3, 587.56e-9, num_rays=500, num_rings=10)
solve_system!(system, source)

# Spot-Diagramm auswerten
hits = spot_diagram(detector)
println("Treffer auf dem Detektor: ", length(hits))
```

### 2. Standalone Julia-Code generieren

```julia
using BeamletOpticsZMX

# Erzeugt ein lesbares Julia-Skript
generate_bmo_script("input_lens.zmx", "reconstructed_lens.jl")
```

Das generierte Skript sieht beispielsweise so aus:

```julia
using BeamletOptics
using LinearAlgebra

const mm = 1e-3

# Lens_1 (Glass: N-BAK1, n ≈ 1.5725)
lens_1 = SphericalLens(14.16mm, 69.41mm, 1.386mm, 12.6mm, λ -> 1.5725)
translate3d!(lens_1, [0.0, 4.5mm, 0.0])

# Doublet_2 (Cemented: N-BAK1 + N-BALF4)
doublet_2 = SphericalDoubletLens(16.2mm, -19.66mm, 8.793mm, 2.313mm, 1.386mm, 11.88mm, λ -> 1.5725, λ -> 1.57956)
translate3d!(doublet_2, [0.0, 6.003mm, 0.0])

# --- Detector ---
detector = Detector(45.23mm)
translate3d!(detector, [0.0, 61.88mm, 0.0])

optics = ObjectGroup([lens_1, doublet_2, ...])
system = System([optics, detector])
```

### 3. CLI Tool

```bash
# System-Metadaten und Baugruppen anzeigen:
julia --project=. bin/zmx2bmo.jl path/to/lens.zmx --info

# ZMX in Julia-Skript konvertieren:
julia --project=. bin/zmx2bmo.jl path/to/lens.zmx output_lens.jl
```

---

## Test-Suite

Die Test-Suite verifiziert alle Teilsysteme sowie End-to-End-Raytracing und Chargen-Tests gegen 35 reale Zemax-Dateien aus Dan Reileys Patent-Bibliothek:

```bash
julia --project=. test/runtests.jl
```

