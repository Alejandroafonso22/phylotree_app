import { HttpClient } from '@angular/common/http';
import { Component, OnInit } from '@angular/core';


export class Species {
  constructor(
    public specie_id: number,
    public scientific_name: string,
    public colloquial_name: string,
    public taxon_id: string,
    public image_specie: string,
    public user: number,
  ) {

  }
}

@Component({
  selector: 'app-species',
  templateUrl: './species.component.html',
  styleUrls: ['./species.component.css']
})
export class SpeciesComponent implements OnInit {
  species: Species[];
  constructor(
    private httpClient: HttpClient) {

  }

  ngOnInit(): void {
    this.getSpecies
  }
  getSpecies() {
    this.httpClient.get<any>('http://192.168.1.33:8000/api/species/?format=json').subscribe(
      response => {
        console.log(response)
        this.species = response;
      }
    )
  }

}
