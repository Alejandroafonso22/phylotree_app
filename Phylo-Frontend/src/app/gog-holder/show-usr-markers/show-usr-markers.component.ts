import { Component, OnInit, Input, OnChanges, SimpleChanges } from '@angular/core';

// Services
import { MarkersService } from '../../services/gog/markers/markers.service';

// Models
import { Marker } from '../../models/marker';

@Component({
  selector: 'app-show-usr-markers',
  templateUrl: './show-usr-markers.component.html',
  styleUrls: ['./show-usr-markers.component.css']
})
export class ShowUsrMarkersComponent implements OnInit, OnChanges {

  // List of markers
  @Input() markersList: Marker[];
  @Input() markerQueryStatus: boolean;
  @Input() species: Array<any>;

  // List of markers that accomplish filter criterias
  filteredMarkersList: Marker[] = [];

  // Current page of pagination
  public currentPage: number;

  // Number of elements in each page
  public itemsPerPage: number = 5;

  // Number of markers
  public markersNumber: number;

  // Filter properties
  public countryCodeFilter: string;
  public dateBeforeFilter: Date;
  public dateAfterFilter: Date;
  public sciNameFilter: string;
  public collNameFilter: string;

  constructor(private service: MarkersService) { }

  ngOnInit(): void {
    this.setSpeciesNames();
    this.refreshMarkersList();
  }

  refreshMarkersList() {
      Object.assign(this.filteredMarkersList, this.markersList);
      this.filterMarkers();
  }

  setSpeciesNames() {
    this.markersList.forEach(marker => {
      for (let index = 0; index < this.species.length; index++) {
        if(marker.specie == this.species[index].specie_id) {
          marker.scientific_name = this.species[index].scientific_name;
          marker.colloquial_name = this.species[index].colloquial_name;
        }
      }
    });
  }

  filterMarkers() {
    this.filteredMarkersList = this.markersList.
    filter(marker => {
      let validCountryCode: boolean = false;
      let validDateBefore: boolean = false;
      let validDateAfter: boolean = false;

      if (this.countryCodeFilter && this.countryCodeFilter != "") {
        if (marker.country.toLowerCase().indexOf
          (this.countryCodeFilter.toLowerCase()) != -1) {
            validCountryCode = true;
        }
      } else {
        validCountryCode = true;
      }

      if(this.dateBeforeFilter) {
        if(new Date(marker.date) <= new Date(this.dateBeforeFilter)) {
          validDateBefore = true;
        }
      } else {
        validDateBefore = true;
      }

      if(this.dateAfterFilter) {
        if(new Date(marker.date) >= new Date(this.dateAfterFilter)) {
          validDateAfter = true;
        }
      } else {
        validDateAfter = true;
      }

      return validCountryCode && validDateBefore && validDateAfter;
    })
    this.currentPage = 1;
    this.markersNumber = this.filteredMarkersList.length;
  }

  ngOnChanges(changes: SimpleChanges) {
    if (typeof changes.species !== 'undefined') {
      this.species = changes.species.currentValue;
    }

    if (typeof changes.markerQueryStatus !== 'undefined') {
      this.markerQueryStatus = changes.markerQueryStatus.currentValue;
    }
    
    if (typeof changes.markersList !== 'undefined') {
      this.markersList = changes.markersList.currentValue;
    }
    this.refreshMarkersList();
  }

}
