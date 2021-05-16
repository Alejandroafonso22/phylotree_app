import { Injectable } from '@angular/core';
import { HttpClient, HttpHeaders} from '@angular/common/http';
import { Observable } from 'rxjs';


@Injectable({
  providedIn: 'root'
})
export class DatabaseService {
  readonly APIUrl = "http://192.168.1.33:8000";

  constructor(private http:HttpClient) { }

  getSpecList():Observable<any[]>{
    return this.http.get<any[]>(this.APIUrl + '/api/species/', {
      headers: new HttpHeaders({
        'Authorization': "Token d7daa018ee895e4c78891e9a4f8b59be69435bd1"
        })
    });
  }
}
