import { Component, DoCheck, Input, OnChanges, OnInit, SimpleChanges } from '@angular/core';
import * as L from 'leaflet';

/* Incompatibility with Webpack
  https://github.com/Leaflet/Leaflet/issues/4968
  Aternative solutions:
    - Copy images to Angular assets folder
    - Map node_modules/leaflet/dist/images/ to assets
*/
import "leaflet/dist/images/marker-shadow.png";
import "../../../assets/SliderControl.js";
// Models
import { Marker } from '../../models/marker';
import { Specie } from '../../models/specie';

@Component({
  selector: 'app-map-view',
  templateUrl: './map-view.component.html',
  styleUrls: ['./map-view.component.css']
})
export class MapViewComponent implements OnInit, OnChanges, DoCheck {

  @Input() markersList: Marker[];
  @Input() species: Specie[];

  private finishFlag: boolean = false;

  private map;
  sliderControl;

  constructor() { }

  ngOnInit(): void {
    this.initMap();
  }

  private initMap() {
    this.map = L.map('map').setView([51.966, 7.6], 10);

    L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}', 
    {attribution: 'Source: Esri, i-cubed, USDA, USGS, AEX,GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community',ext:'jpg',
    }).
    addTo(this.map);

  }

  ngOnChanges(changes: SimpleChanges) {
    if (typeof changes.species !== 'undefined') {
      this.species = changes.species.currentValue;
    }

  }

  ngDoCheck() {
    if(this.species.length != 0) {
      let checker = arr => arr.every(v => v.finished == true);

      if(checker(this.species) && !this.finishFlag) {
        this.finishFlag = true;
        let myRenderer = L.canvas({ padding: 0.5 });

        /*for (var i = 0; i < this.markersList.length; i += 1) { // 100k points
          L.circleMarker([this.markersList[i].latitude, this.markersList[i].longitude], {
            renderer: myRenderer
          }).addTo(this.map).bindPopup('marker ' + i);
        }*/

        this.generateGeoJson(this.markersList)
      }
    }
  }

  generateGeoJson(markerArray: Marker[]) {
    let geoJson = {"type":"FeatureCollection", "features": []}
    for (let index = 0; index < markerArray.length; index++) {
      geoJson.features.push({
		"type": "Feature",
		"properties": {
			"time": `${markerArray[index].date}`,
			"popupContent": "This is a B-Cycle Station. Come pick up a bike and pay by the hour. What a deal!"
		},
		"geometry": {
			"type": "Point",
			"coordinates": [markerArray[index].longitude, markerArray[index].latitude, 1]
		}
	})
    }

    this.addMarkers(geoJson);
  }

  addMarkers(geoJson) {

    let myRenderer = L.canvas({ padding: 0.5 })
    for (var i = 0; i < this.markersList.length; i += 1) { // 100k points
      L.circleMarker([this.markersList[i].latitude, this.markersList[i].longitude], {
        renderer: myRenderer
      }).addTo(this.map).bindPopup('marker marker marker marker marker marker marker marker marker' + i);
    }

    var testlayer = L.geoJson(geoJson, {

      style: function (feature) {
        return feature.properties && feature.properties.style;
      },
      
      onEachFeature: this.onEachFeature,
      
      pointToLayer: function (feature, latlng) {
        return L.circleMarker(latlng, {
          radius: 8,
          fillColor: "rgba(255, 255, 255, 0)",
          color: "rgba(255, 0, 0, 1)",
          weight: 2,
          opacity: 0.5,
          fillOpacity: 0.8
        });
      }
      });
        this.sliderControl = L.control.sliderControl({
        position: "topright",
        layer: testlayer,
        range: true
        });
        //Make sure to add the slider to the map ;-)
        this.map.addControl(this.sliderControl);
        //An initialize the slider
        this.sliderControl.startSlider();
    
        let baseMaps = {
          "All time": myRenderer,
          "Time slider": testlayer
        }
        this.map = L.layerGroup([testlayer, myRenderer])
  }

  onEachFeature(feature, layer) {
		var popupContent = "<p>I started out as a GeoJSON " +
				feature.geometry.type + ", but now I'm a Leaflet vector!</p>";

		if (feature.properties && feature.properties.popupContent) {
			popupContent += feature.properties.popupContent;
		}

		layer.bindPopup(popupContent);

	}
}