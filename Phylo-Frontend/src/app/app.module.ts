import { NgModule } from '@angular/core';
import { BrowserModule } from '@angular/platform-browser';
import { AppRoutingModule } from './app-routing.module';
import { AppComponent } from './app.component';
import { FormsModule } from '@angular/forms';

import {HttpClientModule} from '@angular/common/http'
 
// General imports
import { HomeComponent } from './layout/home/home.component';
import { FooterComponent } from './layout/footer/footer.component';
import { NavbarComponent } from './layout/navbar/navbar.component';
import { HeaderComponent } from './layout/header/header.component';
import { NotFoundComponent } from './meta/not-found/not-found.component';

// GOF app imports
import { GofHolderComponent } from './gof-holder/gof-holder.component';
import { MapViewComponent } from './gof-holder/map-view/map-view.component';
import { SpeciesComponent } from './species/species.component';


@NgModule({
  declarations: [
    AppComponent,
    GofHolderComponent,
    HomeComponent,
    FooterComponent,
    NavbarComponent,
    HeaderComponent,
    NotFoundComponent,
    SpeciesComponent,
  ],
  imports: [
    BrowserModule,
    AppRoutingModule,
    FormsModule,
    HttpClientModule
  ],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule { }
