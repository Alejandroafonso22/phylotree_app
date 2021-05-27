import { Injectable } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
import { Observable } from 'rxjs';

@Injectable({
  providedIn: 'root',
})
export class BackendConnectService {
  readonly APIUrl = "http://127.0.0.1:8000";
  readonly ApiurlSeq = "http://127.0.0.1:8000/api/UploadSeqFromFE/";
  readonly ApiurlGenBank = 'http://127.0.0.1:8000/api/PhylotreeApp/Genbank';
  constructor(private http: HttpClient) { }

  getSpecList(): Observable<any[]> {
    return this.http.get<any[]>(this.APIUrl + '/api/species/', {
      headers: new HttpHeaders({
        'Authorization': "Token 2e46bc7651fe1c925c2513db0eec07377f18518d"
      })
    });
  }

  sendSequenceToBackEnd(acc_number, gene, sequence) {
    const url = this.ApiurlSeq + acc_number + "/" + gene + "/" + sequence;
    return this.http.get<any[]>(url, {
      headers: new HttpHeaders({
        'Authorization': "Token 2e46bc7651fe1c925c2513db0eec07377f18518d"
      })
    });
  }
  
  getGenBanksFromBackEnd() {
    return this.http.get<any[]>(this.ApiurlGenBank, {
      headers: new HttpHeaders({
        'Authorization': "Token 2e46bc7651fe1c925c2513db0eec07377f18518d"
      })
    })
  }


}


