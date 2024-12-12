// ------------------------------------------------------------------------------------------------
// Programm "Labelling"
//   Das Programm dient zur Demonstration verschiedener Objektmerkmale
//   - Der Name des Eingangsbildes wird (ohne ".bmp"-Endung)
//     dem Programm als Argument uebergeben.
//     (Bei Code::Blocks werden die Argumente unter
//      "Project"->"Set programs' arguments..." angegeben)
//
// B. Lang, HS Osnabrueck
// Version Dezember 2016 (RH)
// ------------------------------------------------------------------------------------------------

#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <cstring>
#include <cmath>
#include "BmpRead.h"
#include "BmpWrite.h"
#include "ConvImg.h"
#include "Img.h"

#define INVERTED true

#ifndef M_PI
constexpr double M_PI = 4 * atan(-1);
#endif

// ---------------------------------------------------------------------------
// Aus Uebung 4 uebernehmen
// Optimale Schwelle
// ---------------------------------------------------------------------------
// Parameter:
// [in]  gray_image : Grauwertbild
// Return:
// Binaerbild
// ---------------------------------------------------------------------------
Img<bool> optimal_threshold(const Img<unsigned char> &gray_image)
{
    const unsigned int Width = gray_image.Width();
    const unsigned int Height  = gray_image.Height();
    Img<bool> binary_image(Width, Height);

    // Histogramm berechnen
    std::vector<unsigned int> histogram(256, 0);
    for (unsigned int y = 0; y < Height; ++y) {
        for (unsigned int x = 0; x < Width; ++x) {
            ++histogram[gray_image[y][x]];
        }
    }

    // Gesamtanzahl der Pixel
    const unsigned int total_pixels = Width * Height;

    // Otsu-Methode: Optimalen Schwellenwert berechnen
    double sum = 0.0;
    for (unsigned int i = 0; i < 256; ++i) {
        sum += i * histogram[i];
    }

    double sum_background = 0.0;
    unsigned int weight_background = 0;
    unsigned int weight_foreground = 0;

    double max_variance = 0.0;
    unsigned int optimal_threshold = 0;

    for (unsigned int t = 0; t < 256; ++t) {
        weight_background += histogram[t];
        if (weight_background == 0) continue;

        weight_foreground = total_pixels - weight_background;
        if (weight_foreground == 0) break;

        sum_background += t * histogram[t];

        double mean_background = sum_background / weight_background;
        double mean_foreground = (sum - sum_background) / weight_foreground;

        double variance_between = weight_background * weight_foreground *
                                  std::pow(mean_background - mean_foreground, 2);

        if (variance_between > max_variance) {
            max_variance = variance_between;
            optimal_threshold = t;
        }
    }

    // Binärbild erzeugen
    for (unsigned int y = 0; y < Height; ++y) {
        for (unsigned int x = 0; x < Width; ++x) {
            binary_image[y][x] = gray_image[y][x] > optimal_threshold;
        }
    }

    return binary_image;
}

// ---------------------------------------------------------------------------
// Aus Uebung 4 uebernehmen
// Vektor mit relativen Positionen der Pixel eines quadratischen SEs fuellen
// ---------------------------------------------------------------------------
// Parameter:
// [in]  Diameter : Seitenlaenge des quadratischen SEs in Pixel
// Return:
// Vektor mit relativen Positionen des SEs
// ---------------------------------------------------------------------------
vector<Position> create_square_SE(int Diameter)
{
    vector<Position> ImageWindow;
    int Radius = Diameter / 2;
    for (int dy = -Radius; dy <= Radius; dy++) {
        for (int dx = -Radius; dx <= Radius; dx++) {
            ImageWindow.push_back(Position(dx, dy));
        }
    }
    return ImageWindow;
}

// ---------------------------------------------------------------------------
// Aus Uebung 4 uebernehmen
// Erosion
// ---------------------------------------------------------------------------
// Parameter:
// [in]  src         : Eingangsbild
// [in]  ImageWindow : Strukturierendes Element (Vektor mit Positionen)
// Return:
// erodiertes Bild
// ---------------------------------------------------------------------------
template<typename Pixel>
Img<Pixel> erode(const Img<Pixel> &src, const vector<Position> &ImageWindow)
{
    Img<Pixel> eroded(src.Width(), src.Height());
    const unsigned int Width = src.Width();
    const unsigned int Height = src.Height();


    for (unsigned int y = 0; y < Height; ++y) {
        for (unsigned int x = 0; x < Width; ++x) {
            Pixel min_value = std::numeric_limits<Pixel>::max();
            for (const auto &pos : ImageWindow) {
                int nx = x + pos.get_x();
                int ny = y + pos.get_y();

                if (nx >= 0 && nx < Width && ny >= 0 && ny < Height) {
                    min_value = std::min(min_value, src[ny][nx]);
                }
            }
            eroded[y][x] = min_value;
        }
    }

    return eroded;
}
// ---------------------------------------------------------------------------
// Aus Uebung 4 uebernehmen
// Dilation
// ---------------------------------------------------------------------------
// Parameter:
// [in]  src         : Eingangsbild
// [in]  ImageWindow : Strukturierendes Element (Vektor mit Positionen)
// Return:
// dilatiertes Bild
// ---------------------------------------------------------------------------
template<typename Pixel>
Img<Pixel> dilate(const Img<Pixel> &src, const vector<Position> &ImageWindow)
{
    const unsigned int Width = src.Width();
    const unsigned int Height = src.Height();
    Img<Pixel> dilated(src.Width(), src.Height());


    for (unsigned int y = 0; y < Height; ++y) {
        for (unsigned int x = 0; x < Width; ++x) {
            Pixel max_value = std::numeric_limits<Pixel>::min();
            for (const auto &pos : ImageWindow) {
                int nx = x + pos.get_x();
                int ny = y + pos.get_y();

                if (nx >= 0 && nx < Width && ny >= 0 && ny < Height) {
                    max_value = std::max(max_value, src[ny][nx]);
                }
            }
            dilated[y][x] = max_value;
        }
    }

    return dilated;
}
// ---------------------------------------------------------------------------
// Aus Uebung 4 uebernehmen
// Bilddifferenz
// ---------------------------------------------------------------------------
// Parameter:
// [in]  l : Bild links vom Operator '-'
// [in]  r : Bild rechts vom Operator '-'
// Return:
// Differenzbild
// ---------------------------------------------------------------------------
template<typename Pixel>
Img<Pixel> operator-(const Img<Pixel> &l, const Img<Pixel> &r)
{
    Img<Pixel> d(l.Width(), l.Height());


    for (unsigned int y = 0; y < l.Height(); ++y) {
        for (unsigned int x = 0; x < l.Width(); ++x) {
            // Berechnung der Differenz der Pixelwerte (Clamping bei Bedarf)
            d[y][x] = std::max(0, l[y][x] - r[y][x]);
        }
    }

    return d;
}

// ---------------------------------------------------------------------------
// Aufgabe 2:
// Anzahl der Randpixel pro Objekt berechnen
// Objekte mit dem Label 0 (Hintergrund) werden nicht bearbeitet.
// Im Ausgangsvektor findet sich unter Index N das Ergebnis fuer das
// Objekt mit dem Label N+1.
// ---------------------------------------------------------------------------
// Parameter:
// [in]  label_image : Bild mit Labelnummern (0 = Hintergrund)
// [in]  num_objects : Anzahl der im Bild vorhanden Objekte
// Return:
// Vektor[0..num_objects-1] mit Anzahl der Randpixel je Objekt
// ---------------------------------------------------------------------------
vector<unsigned int> count_MarginPixels(const Img<unsigned int>& label_image, unsigned int num_objects)
{
    unsigned int width = label_image.Width();
    unsigned int height = label_image.Height();

    Img<unsigned int> inner_gradient_label_image;

    // Erosion
    for (unsigned int y = 1; y < height - 1; ++y) {
        for (unsigned int x = 1; x < width - 1; ++x) {
            unsigned int current_label = label_image[y][y];

            if (current_label == 0) continue; // Hintergrund ignorieren

            // Minimalwert im 3x3-Nachbarschaftsfenster berechnen
            unsigned int min_label = current_label;
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dx = -1; dx <= 1; ++dx) {
                    min_label = std::min(min_label, label_image[y+dy][x+dx]);
                }
            }


            inner_gradient_label_image[y][x] = current_label - min_label;
        }
    }


    vector<unsigned int> Pixels(num_objects, 0);

    // Histogramm
    for (unsigned int y = 0; y < height; ++y) {
        for (unsigned int x = 0; x < width; ++x) {
            unsigned int label = inner_gradient_label_image[y][x];
            if (label > 0 && label <= num_objects) {
                Pixels[label - 1]++;
            }
        }
    }

	return Pixels;
}

// ---------------------------------------------------------------------------
// Aufgabe 3:
// Kompaktheit fuer jedes Objekt berechnen
// Objekte mit dem Label 0 (Hintergrund) werden nicht bearbeitet.
// Im Ausgangsvektor findet sich unter Index N das Ergebnis fuer das
// Objekt mit dem Label N+1.
// ---------------------------------------------------------------------------
// Parameter:
// [in]  margin_pixels : Vektor[0..num_objects-1] mit Anzahl der Randpixel je Objekt
// [in]  object_sizes :  Vektor mit Objektflaechen (Anzahl der Pixel) je Objekt
// Return:
// Vektor[0..num_objects-1] mit Kompaktheit je Objekt
// ---------------------------------------------------------------------------
vector<double> compute_Compactness(const vector<unsigned int> &margin_pixels, const vector<unsigned int> &object_sizes)
{
    const unsigned int num_objects = margin_pixels.size();
    vector<double> Compactness(num_objects);

    for (unsigned int i = 0; i < num_objects; ++i) {
        if (object_sizes[i] > 0) {
            Compactness[i] =  (margin_pixels[i] * margin_pixels[i]) / static_cast<double>(object_sizes[i]);
        } else {
            Compactness[i] = 0.0;
        }
    }

    return Compactness;
}

// ------------------------------------------------------------------------------
// Aufgabe 4:
// Exzentrizitaet pro Objekt berechnen
// Objekte mit dem Label 0 (Hintergrund) werden nicht bearbeitet.
// Im Ausgangsvektor findet sich unter Index N das Ergebnis fuer das
// Objekt mit dem Label N+1.
// ---------------------------------------------------------------------------
// Parameter:
// [in]  label_image : Bild mit Labelnummern (0 = Hintergrund)
// [in]  num_objects : Anzahl der im Bild vorhanden Objekte
// Return:
// Vektor[0..num_objects-1] mit Exzentrizitaet je Objekt
// ---------------------------------------------------------------------------
vector<double> compute_Eccentricity(const Img<unsigned int> &LabelImage, unsigned int num_objects)
{
	const unsigned int height = LabelImage.Height();
	const unsigned int width  = LabelImage.Width();
	vector<unsigned int> m00(num_objects);
	vector<unsigned int> m01(num_objects);
	vector<unsigned int> m10(num_objects);
	vector<unsigned int> m11(num_objects);
	vector<unsigned int> m02(num_objects);
	vector<unsigned int> m20(num_objects);
    vector<double> Eccentricities(num_objects, 0.0);

    // Momente berechnen
    for (unsigned int y = 0; y < height; ++y) {
        for (unsigned int x = 0; x < width; ++x) {
            unsigned int label = LabelImage[y][x];
            if (label > 0 && label <= num_objects) {
                unsigned int idx = label - 1;
                m00[idx] += 1;
                m10[idx] += x;
                m01[idx] += y;
                m20[idx] += x * x;
                m02[idx] += y * y;
                m11[idx] += x * y;
            }
        }
    }

    // Exzentrizität für jedes Objekt berechnen
    for (unsigned int i = 0; i < num_objects; ++i) {
        if (m00[i] > 0) {
            double x_mean = static_cast<double>(m10[i]) / m00[i];
            double y_mean = static_cast<double>(m01[i]) / m00[i];

            // Zentrale Momente berechnen
            double mu20 = static_cast<double>(m20[i]) / m00[i] - x_mean * x_mean;
            double mu02 = static_cast<double>(m02[i]) / m00[i] - y_mean * y_mean;
            double mu11 = static_cast<double>(m11[i]) / m00[i] - x_mean * y_mean;

            // Hauptachsen berechnen
            double common = sqrt((mu20 - mu02) * (mu20 - mu02) + 4 * mu11 * mu11);
            double lambda1 = (mu20 + mu02 + common) / 2.0;
            double lambda2 = (mu20 + mu02 - common) / 2.0;

            if (lambda2 > 0) {
                Eccentricities[i] = sqrt(1 - lambda2 / lambda1);
            } else {
                Eccentricities[i] = 0.0;
            }
        } else {
            Eccentricities[i] = 0.0;
        }
    }

    return Eccentricities;
}

// ---------------------------------------------------------------------------
// Labelbild berechnen
// ---------------------------------------------------------------------------
// Parameter:
// [in]  binary_image :
// [in]  r : Bild rechts vom Operator '-'
// Return:
// Differenzbild
// ---------------------------------------------------------------------------
unsigned int Labelling(Img<unsigned int>& label_image, vector <pair<int, int> >& touch_points,
					   vector<unsigned int> &object_sizes, const unsigned int connectivity, const Img<bool>& binary_image)
{
	const Img<bool>& v16(binary_image);
	const unsigned int v0  = v16.Width();
	const unsigned int v1 = v16.Height();
	const unsigned int& v15(connectivity);
	if ((4 != v15) && (8 != v15)) {
		return -1;
	}
	unsigned int v2(0);
	vector <unsigned int> v3;
	Img<unsigned int>&  v12(label_image);
	v12.Resize(v0, v1);
	v12.Margin_Constant(0);
	v3.push_back(0);
	for(unsigned int v4 = 0; v4 < v1; v4++) {
		for(unsigned int v5  = 0; v5 < v0; v5++) {
			v12[v4][v5] = 0;
		}
	}
	vector <pair<int, int> >& v13(touch_points);
	for(unsigned int v4 = 0; v4 < v1; v4++) {
		for(unsigned int v5 = 0; v5 < v0; v5++) {
			if(v16[v4][v5]) {
				vector<unsigned int> v6;
				if(unsigned int v11 = v12[v4 - 1][v5    ]) {
					v6.push_back(v11);
				}
				if(unsigned int v11 = v12[v4    ][v5 - 1]) {
					v6.push_back(v11);
				}
				if(8 == v15) {
					if(unsigned int v11 = v12[v4 - 1][v5 - 1]) {
						v6.push_back(v11);
					}
					if(unsigned int v11 = v12[v4 - 1][v5 + 1]) {
						v6.push_back(v11);
					}
				}
				if(0 == v6.size()) {
					v2++;
					v12[v4][v5] = v2;
					v3.push_back(v2);
					v13.push_back(pair<int, int> (v5, v4));
				} else {
					unsigned int v7 = v6.at(0);
					for(unsigned int v10 = 0; v10 < v6.size(); v10++) {
						unsigned int v8 = v6.at(v10);
						if(v8 < v7) {
							v7 = v8;
						}
					}
					if(v3.at(v7) < v7) {
						v7 = v3.at(v7);
					}
					v12[v4][v5] = v7;
					for(unsigned int v18 = 0; v18 < v6.size(); v18++) {
						unsigned int v8 = v6.at(v18);
						v3.at(v8) = v7;
					}
				}
			}
		}
	}
	for(unsigned int v17 = 0; v17 < v3.size(); v17++) {
		unsigned int v18 = v17;
		while(v18 != v3[v18]) {
			v18 = v3[v18];
		}
		v3[v17] = v18;
	}
	v2 = 0;
	for(unsigned int i = 0; i < v3.size(); i++) {
		if(v3[i] > v2) {
			v2++;
			unsigned int v9 = v3[i];
			for(unsigned int j = i; j < v3.size(); j++) {
				if(v3[j] == v9) {
					v3[j] = v2;
				}
			}
			v13[v2 - 1] = v13[v9 - 1];
		}
	}
	vector<unsigned int>& v14(object_sizes);
	v14.resize(v2, 0);
	for(unsigned int v4 = 0; v4 < v1; v4++) {
		for(unsigned int v5 = 0; v5 < v0; v5++) {
			v12[v4][v5] = v3[v12[v4][v5]];
			if(v12[v4][v5] > 0) {
				v14[v12[v4][v5] - 1] += 1;
			}
		}
	}
	return v2;
}

// ----------------------------------------------------------------------------------------------------
// Labelbild in RGB-Bild konvertieren unter Verwendung eines Vektors,
// der den Labelwerten RGB-Farben zuordnet
// ----------------------------------------------------------------------------------------------------
Img<RGB_Pixel> Labelimage_to_RGB(const Img<unsigned int>& label_image, vector<RGB_Pixel> colors)
{
	const unsigned int width  = label_image.Width();
	const unsigned int height = label_image.Height();
	Img<RGB_Pixel> Label_RGB(width, height);
	for (unsigned int y = 0; y < height; y++) {
		for (unsigned int x = 0; x < width; x++) {
			const unsigned int & label = label_image[y][x];
			if (label == 0) {
				Label_RGB[y][x] = RGB_Pixel(255, 255, 255);
			} else {
				RGB_Pixel &color = colors[label - 1];
				Label_RGB[y][x] = color;
			}
		}
	}
	return Label_RGB;
}

// ----------------------------------------------------------------------------------------------------
// Vektor zufaelliger RGB-Farben fuer die Label eines Labelbilds erzeugen
// ----------------------------------------------------------------------------------------------------
vector<RGB_Pixel> create_LabelColors(unsigned int num_objects)
{
	vector<RGB_Pixel> colors; // enthaelt Farben in Reihenfolge des Farbwinkels

	for (unsigned int i = 0; i < num_objects; i++) {
		double phase((2.0 / 3.0) * 2.0 * M_PI * double(i) / double(num_objects - 1));
		unsigned int blue  = static_cast<unsigned int>(255 * (sin(phase + M_PI / 2.0 + 0.0 * 2 * M_PI / 3.0) + 1.0) / 2.0 + 0.5);
		unsigned int green = static_cast<unsigned int>(255 * (sin(phase + M_PI / 2.0 + 2.0 * 2 * M_PI / 3.0) + 1.0) / 2.0 + 0.5);
		unsigned int red   = static_cast<unsigned int>(255 * (sin(phase + M_PI / 2.0 + 1.0 * 2 * M_PI / 3.0) + 1.0) / 2.0 + 0.5);
		colors.push_back(RGB_Pixel(red, green, blue));
	}

	random_shuffle(colors.begin(), colors.end());
	return colors;
}


// ----------------------------------------------------------------------------------------------------
// Vektor aufsteigender RGB-Farben fuer einen Vektor mit Merkmalswerten erzeugen
// ----------------------------------------------------------------------------------------------------
template<typename VT>
vector<RGB_Pixel> Values2ColorVector(const vector<VT> values)
{
	// Farbe wird abhaengig von den Merkmalswerten vergeben
	vector<VT> sort_values(values);
	sort(sort_values.begin(), sort_values.end());
	VT Max(values[0]);
	for (unsigned int i = 1; i < values.size(); i++) {
		if (values[i] > Max) {
			Max = values[i];
		}
	}
	map<VT, unsigned int> size_to_int;
	vector<RGB_Pixel> sort_colors; // enthaelt Farben in Reihenfolge der Objektgroessen
	for (unsigned int i = 0; i < sort_values.size(); i++) {
		double phase((2.0 / 3.0) * 2.0 * M_PI * double(sort_values[i]) / double(Max));
		unsigned int blue  = static_cast<unsigned int>(255 * (sin(phase + M_PI / 2.0 + 0.0 * 2 * M_PI / 3.0) + 1.0) / 2.0 + 0.5);
		unsigned int green = static_cast<unsigned int>(255 * (sin(phase + M_PI / 2.0 + 2.0 * 2 * M_PI / 3.0) + 1.0) / 2.0 + 0.5);
		unsigned int red   = static_cast<unsigned int>(255 * (sin(phase + M_PI / 2.0 + 1.0 * 2 * M_PI / 3.0) + 1.0) / 2.0 + 0.5);
		sort_colors.push_back(RGB_Pixel(red, green, blue));
		size_to_int[sort_values[i]] = i;
	}
	// Farben werden nun nach Labelnummern sortiert
	vector<RGB_Pixel> colors;
	for (unsigned int i = 0; i < values.size(); i++) {
		colors.push_back(sort_colors[size_to_int[values[i]]]);
	}
	return colors;
}

// ------------------------------------------------------------------------------------------------
// Hauptprogramm
// ------------------------------------------------------------------------------------------------
// Parameter:
// [in] argv[1] : Name des Eingangsbildes
// ------------------------------------------------------------------------------------------------
int main(int argc, char*argv[])
{
	string filename;
	cout << "BV-Praktikum: Labelling" << endl << endl;

	// --------------------------------------------------------------------------------
	// Bildaufnahme
	// --------------------------------------------------------------------------------

	//string Bildname("C:/Users/quint/CLionProjects/labelling/Vorlesungsbeispiel");
	string Bildname("C:/Users/quint/CLionProjects/labelling/Vibrio_cholerae");

	// Bild einlesen
	Img<RGB_Pixel> rgb;
	try {
		filename = Bildname + ".bmp";
		BmpRead(filename.c_str()) >> rgb;
		cout << "Lese " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Lesen von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	// --------------------------------------------------------------------------------
	// Binaeres Quellbild erzeugen
	// --------------------------------------------------------------------------------
	const unsigned int height = rgb.Height();
	const unsigned int width  = rgb.Width();

	Img<unsigned char> uc = ConvImg<unsigned char, RGB_Pixel>(rgb);
	try {
		filename = Bildname + "_uc.bmp";
		BmpWrite(filename.c_str(), uc);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	Img<bool> Binaerbild = optimal_threshold(uc); // Funktion "optimal_threshold" aus Uebung 3 verwenden
	try {
		filename = Bildname + "_bool.bmp";
		BmpWrite(filename.c_str(), Binaerbild);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

#if true == INVERTED
	// Bei Bedarf Binaerbild invertieren: Objektpixel muessen "true" sein
	for (unsigned int y = 0; y < height; y++) {
		for (unsigned int x = 0; x < width; x++) {
			bool &p = Binaerbild[y][x];
			p = not p;
		}
	}
#endif

	// --------------------------------------------------------------------------------
	// Referenzbild mit den Farben zum Einfaerben der Merkmale erzeugen
	// --------------------------------------------------------------------------------
	const unsigned int H(40), B(500);
	vector<unsigned int> ReferenzWerte;
	for (unsigned int i = 0; i < B; i++) {
		ReferenzWerte.push_back(i);
	}
	vector<RGB_Pixel> ReferenzFarben = Values2ColorVector<unsigned int>(ReferenzWerte);
	Img<RGB_Pixel> Referenzbild(B, H);
	for (unsigned int y = 0; y < H; y++) {
		for (unsigned int x = 0; x < B; x++) {
			Referenzbild[y][x] = ReferenzFarben[x];
		}
	}
	try {
		filename = "Referenzbild.bmp";
		BmpWrite(filename.c_str(), Referenzbild);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	// ------------------------------------------------
	// Zu Aufgabe 1: Labelling des Bildes durchfuehren
	// ------------------------------------------------

	Img<unsigned int>       Labelbild;
	vector <pair<int, int> > Antastpunkte;
	vector<unsigned int>    Objektgroessen;
	int Objekte = Labelling(Labelbild, Antastpunkte, Objektgroessen, 8, Binaerbild);
	// Fehlebehandlung
	if (Objekte < 0) {
		cerr << "Fehler beim Labelling" << endl;
		return -1;
	} else if (Objekte == 0) {
		cerr << "Keine Objekte gefunden" << endl;
		return -1;
	}
	unsigned int num_objects = Objektgroessen.size();
	cout << "Gefundene Objekte: " << num_objects << endl;

	{
		// Labelbild mit verschiedenen Farben fuer die Objekte erzeugen und ausgeben
		vector<RGB_Pixel> Farben = create_LabelColors(num_objects);
		Img<RGB_Pixel> LabelAnzeige = Labelimage_to_RGB(Labelbild, Farben);
		for (unsigned int i = 0; i < Objektgroessen.size(); i++) { // Antastpunkte schwarz einzeichnen
			LabelAnzeige[Antastpunkte[i].second][Antastpunkte[i].first] = RGB_Pixel(0, 0, 0);
		}
		try {
			filename = Bildname + "_Labelbild.bmp";
			BmpWrite(filename.c_str(), LabelAnzeige);
			cout << "Schreibe " << filename << endl;
		} catch (const char * s) {
			cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
			return -1;
		}
	}

	// Objektbild erzeugen, in dem die Objekte gemaess ihrer Groesse eingefaerbt sind
	// (Einfaerbung 0:blau, Mitte:gruen, Maximum:rot)
	vector<RGB_Pixel> GroesseFarben = Values2ColorVector<unsigned int>(Objektgroessen);
	Img<RGB_Pixel> GroesseAnzeige = Labelimage_to_RGB(Labelbild, GroesseFarben);
	try {
		filename = Bildname + "_Groesse.bmp";
		BmpWrite(filename.c_str(), GroesseAnzeige);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	// ------------------------------------------------------
	// Zu Aufgabe 2: Anzahl der Randpixel fuer alle Objekte berechnen
	// ------------------------------------------------------

	// Randpixel aller Objekte zaehlen
	vector<unsigned int> RandPixel = count_MarginPixels(Labelbild, num_objects);
	if (RandPixel.size() != Objektgroessen.size()) {
		cerr << "Fehler in der Groesse der Vektoren mit Objektmerkmalen" << endl;
		return -1;
	}
	// Objektbild mit Einfaerbung gemaess Anzahl der Randpixel erzeugen
	vector<RGB_Pixel> RandPixelFarben;
	RandPixelFarben = Values2ColorVector<unsigned int>(RandPixel);
	Img<RGB_Pixel> RandPixelAnzeige = Labelimage_to_RGB(Labelbild, RandPixelFarben);
	try {
		filename = Bildname + "_Randpixel.bmp";
		BmpWrite(filename.c_str(), RandPixelAnzeige);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	// ------------------------------------------------------
	// Zu Aufgabe 3: Kompaktheit fuer alle Objekte berechnen
	// ------------------------------------------------------

	vector<double> Kompaktheiten = compute_Compactness(RandPixel, Objektgroessen);
	vector<RGB_Pixel> KompaktFarben = Values2ColorVector<double>(Kompaktheiten);
	Img<RGB_Pixel> KompaktAnzeige = Labelimage_to_RGB(Labelbild, KompaktFarben);
	try {
		filename = Bildname + "_Kompakt.bmp";
		BmpWrite(filename.c_str(), KompaktAnzeige);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	// ------------------------------------------------------
	// Aufgabe 4: Exzentrizitaet fuer alle Objekte berechnen
	// ------------------------------------------------------
	vector<double> Exzentrizitaeten = compute_Eccentricity(Labelbild, Objektgroessen.size());
	vector<RGB_Pixel> ExzentrizitaetenFarben = Values2ColorVector<double>(Exzentrizitaeten);
	Img<RGB_Pixel> ExzentrizitaetenAnzeige = Labelimage_to_RGB(Labelbild, ExzentrizitaetenFarben);
	try {
		filename = Bildname + "_Exzentrizitaeten.bmp";
		BmpWrite(filename.c_str(), ExzentrizitaetenAnzeige);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	return 0;
}
