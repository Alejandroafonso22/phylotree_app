import { Component, OnInit } from '@angular/core';
import { ShowUsrMarkersComponent } from './show-usr-markers/show-usr-markers.component';
import { MarkersService } from '../services/gog/markers/markers.service';
import { Marker } from '../models/marker';

@Component({
  selector: 'app-gog-holder',
  templateUrl: './gog-holder.component.html',
  styleUrls: ['./gog-holder.component.css']
})
export class GogHolderComponent implements OnInit {

  // List of markers
  markersList: Marker[] = [];

  // Markers observable status
  markerQueryStatus: boolean = false;

  // Error flag
  markerQueryError: boolean = false;

  // Species array
  targetSpecieReady: boolean = false;
  defaultSpeciesReady: boolean = false;
  species: Array<any> = [];

  // User id
  userId: number = 1;
  targetSpecieId: number = 13;

  constructor(private service: MarkersService) { }

  ngOnInit(): void {
    this.cGetTargetSpecie(this.targetSpecieId);
    this.cGetDefaultSpecies();
  }

  cPrepareQueryMarkers() {
    if(this.targetSpecieReady && this.defaultSpeciesReady) {
      for (let index = 0; index < this.species.length; index++) {
        this.species[index].finished =  false;
        this.species[index].empty =  false;
        this.species[index].error =  false;
      }
      
      for (let index = 0; index < this.species.length; index++) {
        this.cGetMarkersList(this.species[index].specie_id, this.userId, index);
      }
    }
  }

  cGetTargetSpecie(specie_id) {
    this.service.getTargetSpecie(specie_id).subscribe(data => {
      this.species = this.species.concat(data);
      this.targetSpecieReady = true;
      this.cPrepareQueryMarkers();
    },
    error => {
      this.markerQueryError = true;
    })
  }

  cGetDefaultSpecies() {
    this.service.getDefaultSpecies().subscribe(data => {
      this.species = this.species.concat(data);
      this.defaultSpeciesReady = true;
      this.cPrepareQueryMarkers();
    },
    error => {
      this.markerQueryError = true;
    })
  }

  cGetMarkersList(specie_id: number, user_id: number, localArrayIndex: number) {
    this.service.getMarkersList(specie_id, user_id).subscribe( data => {
      this.markersList = this.markersList.concat(data);
      if(data.length != 0) {
        this.species[localArrayIndex].empty = false;
      } else {
        this.species[localArrayIndex].empty = true;
      }
      this.species[localArrayIndex].finished = true;
    },
    error => {
      this.species[localArrayIndex].finished = true;
      this.species[localArrayIndex].empty = true;
      this.species[localArrayIndex].error = true;
      this.markerQueryError = true;
    },
    () => {
      this.markerQueryStatus = true;
    })
  }
}
