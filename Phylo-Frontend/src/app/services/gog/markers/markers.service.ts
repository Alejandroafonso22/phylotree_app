import { Injectable } from '@angular/core';

// API consumption imports
import { HttpClient } from '@angular/common/http';
import { Observable } from 'rxjs';
import { catchError } from 'rxjs/operators';

// Models
import { Marker } from '../../../models/marker';

@Injectable({
  providedIn: 'root'
})
export class MarkersService {
  readonly APIUrl = "http://192.168.1.137:8000/api";
  readonly MediUrl = "http://192.168.1.137:8000/api/mapUpload";

  constructor(private http: HttpClient) { }

  getMarkersList(specie_id: number, user_id: number): Observable<Marker[]> {
    return this.http.get<Marker[]>(this.APIUrl + `/markers_with_names/${specie_id}/${user_id}`);
    //.pipe(catchError(this.handleError<any>('getAllFriends', [])));
  }

  getDefaultSpecies(): Observable<any[]> {
    return this.http.get<any[]>(this.APIUrl + `/species/default`);
  }

  getTargetSpecie(specie_id): Observable<any> {
    return this.http.get<any[]>(this.APIUrl + `/species/${specie_id}`)
  }

  addMarker(newMarker: Marker) {
    return this.http.post<any[]>(this.APIUrl + `/markers/`, newMarker);
  }

  updateMarker(targetMarker: Marker) {
    return this.http.put<any[]>(this.APIUrl + `/markers/`, targetMarker);
  }

  deleteMarker(targetMarkerId: number) {
    return this.http.delete<any[]>(this.APIUrl + `/markers/${targetMarkerId}`);
  }

  uploadMapImage(image: any) {
    return this.http.post<any[]>(this.APIUrl + `/markers/save_map`, image);
  }

  private handleError<Marker> (operation: string = "operation", result?: any) {
    return(error: any): Observable<any> => {
      return result;
    }
  }
}
